# This file is part of Bioy
#
#    Bioy is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Bioy is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Bioy.  If not, see <http://www.gnu.org/licenses/>.

"""DEPRECATED: use the classifier subcommand
Classify sequences by grouping blast output by matching taxonomic names

Optional grouping by specimen and query sequences
"""
import sys
import logging

from csv import DictReader, DictWriter
from collections import defaultdict
from math import ceil
from operator import itemgetter

from bioy_pkg import sequtils
from bioy_pkg.utils import Opener, opener, Csv2Dict, groupbyl

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('blast_file',
            nargs = '?',
            default = sys.stdin,
            type = Opener('r'),
            help = 'CSV tabular blast file of query and subject hits')
    parser.add_argument('--all-one-group',
            dest = 'all_one_group',
            action = 'store_true',
            help = """If --map is not provided, the default behavior is to treat
                    all reads as one group; use this option to treat
                    each read as a separate group [%(default)s]""")
    parser.add_argument('-a', '--asterisk',
            default = 100,
            metavar='PERCENT',
            type = float,
            help = 'Next to any species above a certain threshold [%(default)s]')
    parser.add_argument('--copy-numbers',
            metavar = 'CSV',
            type = Opener(),
            help = 'columns: tax_id, median')
    parser.add_argument('-c', '--coverage',
            default = 95,
            metavar = 'PERCENT',
            type = float,
            help = 'percent of alignment coverage of blast result [%(default)s]')
    parser.add_argument('--details-identity',
            metavar = 'PERCENT',
            help = 'Minimum identity to include blast hits in details file',
            type = float,
            default = 90)
    parser.add_argument('--details-full',
            action = 'store_true',
            help = 'do not limit out_details to only larget cluster per assignment')
    parser.add_argument('--exclude-by-taxid',
            metavar = 'CSV',
            type = lambda f: set(e for e in DictReader(opener(f), fieldnames ='tax_id')),
            default = {},
            help = 'column: tax_id')
    parser.add_argument('--group-def',
            metavar = 'INT',
            action = 'append',
            default = [],
            help = """define a group threshold for a particular rank overriding
                      --target-max-group-size. example: genus:2""")
    parser.add_argument('--group-label',
            metavar = 'LABEL',
            default = 'all',
            help = 'Single group label for reads')
    parser.add_argument('-o', '--out',
            default = sys.stdout,
            type = Opener('w'),
            metavar = 'CSV',
            help = """columns: specimen, max_percent, min_percent, max_coverage,
                      min_coverage, assignment_id, assignment, clusters, reads,
                      pct_reads, corrected, pct_corrected, target_rank, hi, low, tax_ids""")
    parser.add_argument('-m', '--map',
            metavar = 'CSV',
            type = Opener(),
            default = {},
            help = 'columns: name, specimen')
    parser.add_argument('--max-ambiguous',
            metavar = 'INT',
            default = 3,
            type = int,
            help = 'Maximum ambiguous count in reference sequences [%(default)s]')
    parser.add_argument('--max-identity',
            default = 100,
            metavar = 'PERCENT',
            type = float,
            help = 'maximum identity threshold for accepting matches [<= %(default)s]')
    parser.add_argument('--min-cluster-size',
            default = 0,
            metavar = 'INT',
            type = int,
            help = 'minimum cluster size to include in classification output')
    parser.add_argument('--min-identity',
            default = 99,
            metavar = 'PERCENT',
            type = float,
            help = 'minimum identity threshold for accepting matches [> %(default)s]')
    parser.add_argument('-s', '--seq-info',
            required = True,
            metavar = 'CSV',
            type = Opener(),
            help = 'seq info file(s) to match sequence ids to taxids [%(default)s]')
    parser.add_argument('-t', '--taxonomy',
            required = True,
            metavar = 'CSV',
            type = Csv2Dict('tax_id'),
            help = 'tax table of taxids and species names [%(default)s]')
    parser.add_argument('-O', '--out-detail',
            type  = lambda f: DictWriter(opener(f, 'w'), extrasaction = 'ignore', fieldnames = [
                'specimen', 'assignment', 'assignment_id', 'qseqid', 'sseqid', 'pident', 'coverage', 'ambig_count',
                'accession', 'tax_id', 'tax_name', 'target_rank', 'rank', 'hi', 'low'
                ]),
            metavar = 'CSV',
            help = """columns: specimen, assignment, assignment_id,
                      qseqid, sseqid, pident, coverage, ambig_count,
                      accession, tax_id, tax_name, target_rank, rank, hi, low""")
    parser.add_argument('--target-max-group-size',
            metavar = 'INTEGER',
            default = 3,
            type = int,
            help = """group multiple target-rank assignments that
                      excede a threshold to a higher rank [%(default)s]""")
    parser.add_argument('--target-rank',
            metavar='RANK',
            help = 'Rank at which to classify. Default: "%(default)s"',
            default = 'species')
    parser.add_argument('-w', '--weights',
            metavar = 'CSV',
            type = Opener(),
            help = 'columns: name, weight')
    ### csv.Sniffer.has_header is *not* reliable enough
    parser.add_argument('--has-header', action = 'store_true',
            help = 'specify this if blast data has a header')

def coverage(start, end, length):
    return (float(end) - float(start) + 1) / float(length) * 100

def mean(l):
    l = list(l)
    return float(sum(l)) / len(l) if len(l) > 0 else 0

def condense(queries, floor_rank, max_size, ranks, rank_thresholds, target_rank = None):
    target_rank = target_rank or ranks[0]

    groups = list(groupbyl(queries, key = itemgetter(target_rank)))

    num_groups = len(groups)

    if rank_thresholds.get(target_rank, max_size) < num_groups:
        return queries

    # assign where available target_rank_ids
    # groups without 'i' values remain assigned at previous (higher) rank
    for g in (g for i,g in groups if i):
        for q in g:
            q['target_rank_id'] = q[target_rank]

    # return if we hit the floor
    if target_rank == floor_rank:
        return queries

    # else move down a rank
    target_rank = ranks[ranks.index(target_rank) + 1]

    # recurse down the tax tree
    condensed = []
    for _,g in groups:
        c = condense(g, floor_rank, max_size, ranks, rank_thresholds, target_rank)
        condensed.extend(c)

    return condensed

def action(args):
    ### format format blast data and add additional available information
    fieldnames = None if args.has_header else sequtils.BLAST_HEADER_DEFAULT
    blast_results = DictReader(args.blast_file, fieldnames = fieldnames)
    blast_results = list(blast_results)

    sseqids = set(s['sseqid'] for s in blast_results)
    qseqids = set(s['qseqid'] for s in blast_results)

    # load seq_info and map file
    mapfile = DictReader(args.map, fieldnames = ['name', 'specimen'])
    mapfile = {m['name']:m['specimen'] for m in mapfile if m['name'] in qseqids}

    seq_info = DictReader(args.seq_info)
    seq_info = {s['seqname']:s for s in seq_info if s['seqname'] in sseqids}

    # pident
    def pident(b):
        return dict(b, pident = float(b['pident'])) if b['sseqid'] else b

    blast_results = (pident(b) for b in blast_results)

    # coverage
    def cov(b):
        if b['sseqid'] and b['qcovs']:
            b['coverage'] = float(b['qcovs'])
            return b
        elif b['sseqid']:
            c = coverage(b['qstart'], b['qend'], b['qlen'])
            return dict(b, coverage = c)
        else:
            return b

    blast_results = (cov(b) for b in blast_results)

    # seq info
    def info(b):
        return dict(seq_info[b['sseqid']], **b) if b['sseqid'] else b

    blast_results = (info(b) for b in blast_results)

    # tax info
    def tax_info(b):
        return dict(args.taxonomy[b['tax_id']], **b) if b['sseqid'] else b

    blast_results = (tax_info(b) for b in blast_results)

    ### output file headers
    fieldnames = ['specimen', 'max_percent', 'min_percent', 'max_coverage',
                  'min_coverage', 'assignment_id', 'assignment']

    if args.weights:
        weights = DictReader(args.weights, fieldnames = ['name', 'weight'])
        weights = {d['name']:d['weight'] for d in weights if d['name'] in qseqids}
        fieldnames += ['clusters', 'reads', 'pct_reads']
    else:
        weights = {}

    if args.copy_numbers:
        copy_numbers = DictReader(args.copy_numbers)
        copy_numbers = {d['tax_id']:float(d['median']) for d in copy_numbers}
        fieldnames += ['corrected', 'pct_corrected']
    else:
        copy_numbers = {}

    # TODO: take out target_rank, hi, low and provide in pipeline using csvmod
    # TODO: option to include tax_ids (default no)
    fieldnames += ['target_rank', 'hi', 'low', 'tax_ids']

    ### Columns
    out = DictWriter(args.out,
            extrasaction = 'ignore',
            fieldnames = fieldnames)
    out.writeheader()

    if args.out_detail:
        args.out_detail.writeheader()

    def blast_hit(hit, args):
        return hit['sseqid'] and \
               hit[args.target_rank] and \
               hit['coverage'] >= args.coverage and \
               float(weights.get(hit['qseqid'], 1)) >= args.min_cluster_size and \
               hit[args.target_rank] not in args.exclude_by_taxid and \
               hit['qseqid'] != hit['sseqid'] and \
               int(hit['ambig_count']) <= args.max_ambiguous

    ### Rows
    etc = '[no blast result]' # This row will hold all unmatched

    # groups have list position prioritization
    groups = [
        ('> {}%'.format(args.max_identity),
            lambda h: blast_hit(h, args) and h['pident'] > args.max_identity),
        (None,
            lambda h: blast_hit(h, args) and args.max_identity >= h['pident'] > args.min_identity),
        ('<= {}%'.format(args.min_identity),
            lambda h: blast_hit(h, args) and h['pident'] <= args.min_identity),
    ]

    # used later for results output
    group_cats = map(itemgetter(0), groups)
    group_cats.append(etc)

    # assignment rank thresholds
    rank_thresholds = (d.split(':') for d in args.group_def)
    rank_thresholds = dict((k, int(v)) for k,v in rank_thresholds)

    # rt = {k: int(v) for k, v in (d.split(':') for d in args.group_def)}

    # group by specimen
    if args.map:
        specimen_grouper = lambda s: mapfile[s['qseqid']]
    elif args.all_one_group:
        specimen_grouper = lambda s: args.group_label
    else:
        specimen_grouper = lambda s: s['qseqid']

    blast_results = groupbyl(blast_results, key = specimen_grouper)

    assignments = [] # assignment list for assignment ids

    for specimen, hits in blast_results:
        categories = defaultdict(list)

        # clusters will hold the query ids as hits are matched to categories
        clusters = set()

        # filter out categories
        for cat, fltr in groups:
            matches = filter(fltr, hits)

            if cat:
                categories[cat] = matches
            else:
                # create sets of tax_rank_id
                query_group = groupbyl(matches, key = itemgetter('qseqid'))

                target_cats = defaultdict(list)
                for _,queries in query_group:
                    queries = condense(
                            queries,
                            args.target_rank,
                            args.target_max_group_size,
                            sequtils.RANKS,
                            rank_thresholds)
                    cat = map(itemgetter('target_rank_id'), queries)
                    cat = frozenset(cat)

                    target_cats[cat].extend(queries)

                categories = dict(categories, **target_cats)

            # add query ids that were matched to a filter
            clusters |= set(map(itemgetter('qseqid'), matches))

            # remove all hits corresponding to a matched query id (cluster)
            hits = filter(lambda h: h['qseqid'] not in clusters, hits)

        # remaining hits go in the etc ('no match') category
        categories[etc] = hits

        # calculate read counts
        read_counts = dict()
        for k,v in categories.items():
            qseqids = set(map(itemgetter('qseqid'), v))
            weight = sum(float(weights.get(q, 1)) for q in qseqids)
            read_counts[k] = weight

        taxids = set()
        for k,v in categories.items():
            if k is not etc:
                for h in v:
                    taxids.add(h['tax_id'])

        ### list of assigned ids for count corrections
        assigned_ids = dict()
        for k,v in categories.items():
            if k is not etc and v:
                assigned_ids[k] = set(map(itemgetter('tax_id'), v))

        # correction counts
        corrected_counts = dict()
        for k,v in categories.items():
            if k is not etc and v:
                av = mean(copy_numbers.get(t, 1) for t in assigned_ids[k])
                corrected_counts[k] = ceil(read_counts[k] / av)

        # finally take the root value for the etc category
        corrected_counts[etc] = ceil(read_counts[etc] / copy_numbers.get('1', 1))

        # totals for percent calculations later
        total_reads = sum(v for v in read_counts.values())
        total_corrected = sum(v for v in corrected_counts.values())

        # Print classifications per specimen sorted by # of reads in reverse (descending) order

        sort_by_reads_assign = lambda (c,h): corrected_counts.get(c, None)

        for cat, hits in sorted(categories.items(), key = sort_by_reads_assign, reverse = True):

            # continue if their are hits
            if hits:

                # for incrementing assignment id's
                if cat not in assignments:
                    assignments.append(cat)

                assignment_id = assignments.index(cat)

                reads = read_counts[cat]
                reads_corrected = corrected_counts[cat]

                clusters = set(map(itemgetter('qseqid'), hits))

                results = dict(
                        hi = args.max_identity,
                        low = args.min_identity,
                        target_rank = args.target_rank,
                        specimen = specimen,
                        assignment_id = assignment_id,
                        reads = int(reads),
                        pct_reads = '{0:.2f}'.format(reads / total_reads * 100),
                        corrected = int(reads_corrected),
                        pct_corrected = '{0:.2f}'.format(reads_corrected / total_corrected * 100),
                        clusters = len(clusters))

                if cat is etc:
                    assignment = etc
                    results = dict(results, assignment = assignment)
                else:
                    taxids = set(map(itemgetter('tax_id'), hits))
                    coverages = set(map(itemgetter('coverage'), hits))
                    percents = set(map(itemgetter('pident'), hits))

                    if cat in group_cats:
                        assignment = cat
                    else:
                        names = [args.taxonomy[h['target_rank_id']]['tax_name'] for h in  hits]
                        selectors = [h['pident'] >= args.asterisk for h in hits]
                        assignment = sequtils.format_taxonomy(names, selectors, '*')

                    results = dict(results,
                        assignment = assignment,
                        max_percent = '{0:.2f}'.format(max(percents)),
                        min_percent = '{0:.2f}'.format(min(percents)),
                        max_coverage = '{0:.2f}'.format(max(coverages)),
                        min_coverage = '{0:.2f}'.format(min(coverages)),
                        tax_ids = ' '.join(taxids))

                out.writerow(results)

                if args.out_detail:
                    if not args.details_full:
                        # drop the no_hits
                        hits = [h for h in hits if 'tax_id' in h]
                        # only report heaviest centroid
                        clusters_and_sizes = [(float(weights.get(c, 1.0)), c) for c in clusters]
                        _, largest = max(clusters_and_sizes)
                        hits = (h for h in hits if h['qseqid'] == largest)

                    for h in hits:
                        args.out_detail.writerow(dict(
                        specimen = specimen,
                        assignment = assignment,
                        assignment_id = assignment_id,
                        hi = args.max_identity,
                        low = args.min_identity,
                        target_rank = args.target_rank,
                        **h))

