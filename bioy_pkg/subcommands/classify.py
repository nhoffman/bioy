"""
Classify sequences by grouping blast output by matching taxonomic names

Optional grouping by specimen and query sequences
"""
import sys
import logging

from csv import DictReader, DictWriter
from collections import defaultdict
from itertools import groupby, imap, ifilter
from operator import itemgetter

from bioy_pkg.sequtils import UNCLASSIFIED_REGEX, format_taxonomy
from bioy_pkg.utils import Opener, Csv2Dict

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('blast_file',
            nargs = '?',
            default = sys.stdin,
            type = Opener('r'),
            help = 'CSV tabular blast file of query and subject hits')
    parser.add_argument('-o', '--out',
            default = sys.stdout,
            type = Opener('w'),
            help = """Tab delimited list of classifications by species
                   [specimen, species, reads, clusters, max_percent, min_percent[]""")
    parser.add_argument('-O', '--out-detail',
            type  = Opener('w'),
            help = 'Add detailed csv file')
    parser.add_argument('-s', '--seq-info',
            required = True,
            type = Csv2Dict('seqname'),
            help = 'seq info file(s) to match sequence ids to taxids [%(default)s]')
    parser.add_argument('-t', '--taxonomy',
            required = True,
            type = Csv2Dict('tax_id'),
            help = 'tax table of taxids and species names [%(default)s]')
    parser.add_argument('--min-identity',
            default = 99,
            type = float,
            help = 'minimum identity threshold for accepting matches [> %(default)s]')
    parser.add_argument('--max-identity',
            default = 100,
            type = float,
            help = 'maximum identity threshold for accepting matches [<= %(default)s]')
    parser.add_argument('-c', '--coverage',
            default = 95,
            type = float,
            help = 'percent of alignment coverage of blast result [%(default)s]')
    parser.add_argument('-a', '--asterisk',
            default = 100,
            type = float,
            help = 'Next to any species above a certain threshold [[%(default)s]')
    parser.add_argument('--exclude-by-taxid',
            type = Csv2Dict('tax_id'),
            default = {},
            help = """csv file with column "tax_id" providing a set of tax_ids
                      to exclude from the results""")
    parser.add_argument('-w', '--weights',
            type = Csv2Dict('name', 'weight', fieldnames = ['name', 'weight']),
            default = {},
            help = 'Provides a weight (number of reads) associated with each id')
    parser.add_argument('-m', '--map',
            type = Csv2Dict('name', 'specimen', fieldnames = ['name', 'specimen']),
            default = {},
            help = 'map file with sequence ids to specimen names')
    parser.add_argument('--not-all-one-group',
            dest = 'all_one_group',
            action = 'store_false',
            default = True,
            help = """If --map is not provided, the default behavior is to treat
                    all reads as one group; use this option to treat
                    each read as a separate group [%(default)s]""")
    parser.add_argument('--min-cluster-size',
            default = 0,
            type = int,
            help = 'minimum cluster size to include in classification output')
    parser.add_argument('--blast-fmt',
            help = 'blast header (default: qseqid sseqid pident qstart qend qlen)',
            default = 'qseqid sseqid pident qstart qend qlen')
    parser.add_argument('--target-rank',
            help = 'Rank at which to classify. Default: "%(default)s"',
            default = 'species')
    parser.add_argument('--details-identity',
            help = 'Minimum identity threshold for classification details file',
            type = float,
            default = 90)

def update_blast_results(b, seq_info, taxonomy, target_rank):
    info = seq_info[b['sseqid']]
    tax = taxonomy[info['tax_id']]

    b['query'] = b['qseqid']
    b['subject'] = b['sseqid']
    b['pident'] = float(b['pident'])

    b['tax_id'] = info['tax_id']
    b['accession'] = info['accession']
    b['ambig_count'] = int(info['ambig_count'])
    b['full_name'] = tax['tax_name']
    b['rank'] = tax['rank']

    b['target_rank_id'] = tax[target_rank]

    if b['target_rank_id']:
        b['target_rank_name'] = taxonomy[b['target_rank_id']]['tax_name']

    return b

def coverage(start, end, length):
    return (float(end) - float(start) + 1) / float(length) * 100

def action(args):
    ### Rows
    etc = 'no match' # This row will hold all unmatched

    # groups have list position prioritization
    groups = [
        ('> {}%'.format(args.max_identity), lambda h: h['pident'] > args.max_identity),
        (None, lambda h: args.max_identity >= h['pident'] > args.min_identity),
        ('<= {}%'.format(args.min_identity), lambda h: h['pident'] <= args.min_identity),
    ]

    ### Columns
    out = DictWriter(args.out, extrasaction = 'ignore', fieldnames = [
        'max_percent', 'min_percent', 'max_coverage', 'min_coverage', 'tax_name',
        'specimen', 'reads', 'pct_reads', 'clusters',
        'target_rank', 'hi', 'low',
        ])
    out.writeheader()

    if args.out_detail:
        detail = DictWriter(args.out_detail, extrasaction = 'ignore', fieldnames = [
            'specimen', 'tax_name', 'query', 'subject', 'pident', 'coverage', 'ambig_count',
            'target_rank_id', 'target_rank_name', 'accession', 'tax_id',
            'full_name', 'target_rank', 'hi', 'low'
            ])
        detail.writeheader()
    ###

    ### filter and format format blast data
    blast_fmt = args.blast_fmt.split()

    blast_results = DictReader(args.blast_file, fieldnames = blast_fmt, delimiter = '\t')

    # add some preliminary values to blast results
    blast_results = imap(lambda b: dict({
        'coverage':coverage(b['qstart'], b['qend'], b['qlen']),
        }, **b), blast_results)

    # some raw filtering
    blast_results = ifilter(lambda b:
            float(args.weights.get(b['qseqid'], 1)) >= args.min_cluster_size, blast_results)
    blast_results = ifilter(lambda b: float(b['pident']) >= args.details_identity, blast_results)
    blast_results = ifilter(lambda b: float(b['coverage']) >= args.coverage, blast_results)

    # add required values for classification
    blast_results = imap(lambda b:
            update_blast_results(b, args.seq_info, args.taxonomy, args.target_rank), blast_results)

    # remove hits with no rank ids
    blast_results = ifilter(lambda b: b['target_rank_id'], blast_results)

    blast_results = ifilter(
            lambda b: b['target_rank_id'] not in args.exclude_by_taxid, blast_results)

    # (Optional) In some cases you want to filter some target hits
#    blast_results = ifilter(lambda b: b['qseqid'] != b['sseqid'], blast_results)
#    blast_results = ifilter(lambda b: b['ambig_count'] < 3, blast_results)
#    blast_results = ifilter(
#            lambda b: not UNCLASSIFIED_REGEX.search(b['target_rank_name']), blast_results)
    ###

    # first, group by specimen
    if args.map:
        specimen_grouper = lambda s: args.map[s['qseqid']]
    elif args.all_one_group:
        specimen_grouper = lambda s: 'all'
    else:
        specimen_grouper = lambda s: s['qseqid']

    for specimen, hits in groupby(blast_results, specimen_grouper):
        hits = list(hits)

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
                cats = defaultdict(list)
                for _, queries in groupby(matches, itemgetter('query')):
                    queries = list(queries)
                    target_ids = frozenset(map(itemgetter('target_rank_id'), queries))
                    cats[target_ids].extend(queries)

                # and finally add to categories
                for _, queries in cats.items():
                    queries = list(queries)
                    names = map(itemgetter('target_rank_name'), queries)
                    selectors = map(lambda h: h['pident'] >= args.asterisk, queries)
                    tax = format_taxonomy(names, selectors, '*')
                    categories[tax].extend(queries)

            # add query ids that were matched to a filter
            clusters |= set(map(itemgetter('query'), matches))

            # remove all hits corresponding to a matched query id (cluster)
            hits = filter(lambda h: h['query'] not in clusters, hits)

        # remaining hits go in the etc ('no match') category
        categories[etc] = hits

        # include the remaining hits in the clusters set
        clusters |= set(map(itemgetter('query'), hits))

        total_reads = sum(float(args.weights.get(c, 1)) for c in clusters)

        # Print classifications per specimen sorted by # of reads in reverse (descending) order
        sort_by_reads = lambda c: sum(int(args.weights.get(q, 1)) for q in set(c['query'] for c in c[1]))
        for cat, hits in sorted(categories.items(), key=sort_by_reads, reverse=True):
            # only output categories with hits
            if hits:
                clusters = set(map(itemgetter('query'), hits))
                coverages = set(map(itemgetter('coverage'), hits))
                percents = set(map(itemgetter('pident'), hits))
                reads = sum(float(args.weights.get(c, 1)) for c in clusters)

                out.writerow({
                    'hi':args.max_identity,
                    'low':args.min_identity,
                    'target_rank':args.target_rank,
                    'specimen':specimen,
                    'tax_name':cat,
                    'reads':int(reads),
                    'pct_reads':'{0:.2f}'.format(reads / total_reads * 100),
                    'clusters':len(clusters),
                    'max_percent':'{0:.2f}'.format(max(percents)),
                    'min_percent':'{0:.2f}'.format(min(percents)),
                    'max_coverage':'{0:.2f}'.format(max(coverages)),
                    'min_coverage':'{0:.2f}'.format(min(coverages)),
                    })

                if args.out_detail:
                    detail.writerows(dict({
                        'specimen':specimen,
                        'tax_name':cat,
                        'hi':args.max_identity,
                        'low':args.min_identity,
                        'target_rank':args.target_rank,
                        }, **h) for h in hits)


