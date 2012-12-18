"""
Classify sequences by grouping blast output by matching taxonomic names

Optional grouping by specimen and query sequences
"""
import sys
import logging
import os
import csv
from csv import DictReader, DictWriter
from collections import defaultdict
from itertools import groupby

from ion_tools.sequtils import UNCLASSIFIED_REGEX
from ion_tools.utils import Opener, csv2dict

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('blast_file',
            nargs = '?',
            default = sys.stdin,
            type = Opener('r'),
            help = 'csv tabular blast file of query and subject hits')
    parser.add_argument('-o', '--out',
            default = sys.stdout,
            type = Opener('w'),
            help = 'tab delimited list of classifications by species \
                   [specimen, species, reads, clusters, max_percent, min_percent[]')
    parser.add_argument('-O', '--detail',
            default = os.devnull,
            type  = Opener('w'),
            help = 'add detailed csv file')
    parser.add_argument('-s', '--seq-info',
            required = True,
            type = csv2dict('seqname'),
            help = 'seq info file(s) to match sequence ids to taxids')
    parser.add_argument('-t', '--tax',
            required = True,
            type = csv2dict('tax_id'),
            help = 'tax table of taxids and species names')
    parser.add_argument('-i', '--identity-threshold',
            default = 99,
            type = float,
            help = 'blast percent threshold for classification [%(default)s]')
    parser.add_argument('-c', '--coverage-threshold',
            default = 95,
            type = float,
            help = 'percent of alignment coverage of blast result [%(default)s]')
    parser.add_argument('-a', '--asterisk',
            default = 100,
            type = float,
            help = 'next to any species above a certain threshold')
    parser.add_argument('-F', '--no-filter-by-name',
            action='store_false',
            dest = 'filter_by_name',
            default = True,
            help = 'do not exclude records with unclassified-looking species name')
    parser.add_argument('--exclude-by-taxid',
            type = csv2dict('tax_id'),
            help = 'csv file with column "tax_id" providing a set of tax_ids to exclude from the results')
    parser.add_argument('-w', '--weights',
            type = csv2dict('name', 'weight', ['name', 'weight']),
            default = {},
            help = 'Provides a weight (number of reads) associated with each id')
    parser.add_argument('-m', '--map',
            type = csv2dict('name', 'specimen', ['name', 'specimen']),
            default = {},
            help = 'map file with sequence ids to specimen names')
    parser.add_argument(
            '--not-all-one-group',
            dest = 'all_one_group',
            action = 'store_false',
            default = True,
            help = 'If --map is not provided, the default behavior is to treat \
                    all reads as one group; use this option to treat \
                    each read as a separate group.')
    parser.add_argument(
            '--min-cluster-size',
            default = 0,
            type = int,
            help = 'minimum cluster size to include in classification output')
    parser.add_argument('--no-hits', type = Opener('w'),
                        help = 'file containing name of each query with no hits')
    parser.add_argument('--all-taxids', type = Opener('w'),
                        help = 'file containing set of taxids represented in the output')
    parser.add_argument('--blast-fmt',
            help = 'blast header (default: qseqid sseqid pident qstart qend qlen)',
            default = 'qseqid sseqid pident qstart qend qlen')

def action(args):
    # aggregate(a) blast results and classifications by clusters
    blast_fmt = args.blast_fmt.split()
    out_header = [
        'specimen',
        'reads',
        'pct_reads',
        'clusters',
        'species',
        'max_percent', 'min_percent',
        'max_coverage', 'min_coverage',
        # 'query_ids'
        ]
    detail_header = [
            'specimen',
            'query',
            'subject',
            'match',
            'pident',
            'coverage',
            'ambig',
            'species_name',
            'accession',
            'species_id',
            'tax_id',
            'tax_name',
            'rank'
            ]

    out = DictWriter(args.out, out_header)
    out.writeheader()

    detail = DictWriter(args.detail, detail_header)
    detail.writeheader()

    if args.no_hits:
        args.no_hits.write('query\n')

    blast_results = DictReader(args.blast_file, fieldnames = blast_fmt, delimiter = '\t')

    # first, group by specimen
    if args.map:
        specimen_grouper = lambda s: args.map[s['qseqid']]
    elif args.all_one_group:
        specimen_grouper = lambda s: 'all'
    else:
        specimen_grouper = lambda s: s['qseqid']

    all_taxids = set()
    for specimen, results in groupby(blast_results, specimen_grouper):
        classifications = {}
        # nest, group by query
        s_reads = 0
        for query, hits in groupby(results, lambda q: q['qseqid']):
            q_reads = float(args.weights.get(query, 1))
            if q_reads < args.min_cluster_size:
                continue
            s_reads += q_reads
            data = defaultdict(list)
            for h in hits:
                info = args.seq_info[h['sseqid']]
                taxonomy = args.tax[info['tax_id']]
                species_id = taxonomy['species'] or None
                species_name = args.tax[species_id]['tax_name'] if species_id else None
                coverage = (float(h['qend']) - float(h['qstart']) + 1) / float(h['qlen']) * 100

                match = query != h['sseqid']
                # filter by taxonomy
                if args.filter_by_name:
                    match &= species_id and not UNCLASSIFIED_REGEX.search(species_name)

                match &= query != h['sseqid']
                match &= float(h['pident']) >= args.identity_threshold
                match &= int(info['ambig_count']) < 3
                match &= coverage >= args.coverage_threshold
                if args.exclude_by_taxid:
                    match &= species_id not in args.exclude_by_taxid
                # aggregate data that is a match
                if match:
                    data[species_id].append(dict(h, **{'coverage':coverage}))
                    all_taxids.add((species_id, species_name))

                # Ouput hit details
                if coverage > 90 and float(h['pident']) > 90:
                    detail.writerow({
                        'specimen':specimen,
                        'query':query,
                        'accession': info['accession'],
                        'subject': h['sseqid'],
                        'match':'yes' if match else 'no',
                        'species_name': species_name,
                        'species_id': species_id,
                        'rank': taxonomy['rank'],
                        'tax_id': info['tax_id'],
                        'tax_name': taxonomy['tax_name'],
                        'coverage': coverage,
                        'pident': h['pident'],
                        'ambig': info['ambig_count']
                        })

            # data is empty if there were no hits
            if args.no_hits and not data:
                args.no_hits.write(query+'\n')

            ## Add data into the classification dict
            key = frozenset(data.keys()) if data else frozenset()
            if key in classifications:
                (d, r, c) = classifications[key]
                for k,v in data.items():
                    d[k].extend(v)
                classifications[key] = (d, r + q_reads, c + 1)
            else:
                classifications[key] = (data, q_reads, 1)

        # Print classifications per specimen sorted by # of reads in reverse (descending) order
        for s_ids, (data, reads, clusters) in sorted(
            classifications.items(), key = lambda a: a[1][1], reverse=True):
            # Species names
            if s_ids:
                sp = defaultdict(list)
                for i in s_ids:
                    name = args.tax[i]['tax_name'].split(None, 1)
                    max_pident = max([float(d['pident']) for d in data[i]])
                    max_coverage = max([d['coverage'] for d in data[i]])
                    if max_pident >= args.asterisk and max_coverage >= args.coverage_threshold:
                        star = '*'
                    else:
                        star = ''
                    sp[name[0]].append(star if len(name) == 1 else name[1] + star)
                species = ';'.join(['{} {}'.format(
                    genus,'/'.join(species)) if species else genus for genus,species in sp.items()])
            else:
                species = 'no match'

            percents = [float(v['pident']) for values in data.values() for v in values]
            coverages = [v['coverage'] for values in data.values() for v in values]
            # qseqids = set(v['qseqid'] for values in data.values() for v in values)

            # Output
            out.writerow({
                'specimen':specimen,
                'species':species,
                'reads':int(reads),
                'pct_reads':reads/s_reads * 100,
                'clusters':clusters,
                'max_percent':max(percents) if percents else None,
                'min_percent':min(percents) if percents else None,
                'max_coverage':max(coverages) if coverages else None,
                'min_coverage':min(coverages) if coverages else None,
                # 'query_ids': ' '.join(qseqids)
                })

        if args.all_taxids:
            writer = csv.writer(args.all_taxids)
            writer.writerow(['tax_id','tax_name'])
            writer.writerows(all_taxids)
