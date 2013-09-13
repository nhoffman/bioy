"""
Classify sequences by grouping blast output by matching taxonomic names

Optional grouping by specimen and query sequences
"""
import sys
import logging

from csv import DictReader, DictWriter
from collections import defaultdict
from operator import itemgetter

from bioy_pkg import sequtils
from bioy_pkg.utils import Opener, Csv2Dict, groupbyl

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
            type  = lambda f: DictWriter(Opener('w')(f), extrasaction = 'ignore', fieldnames = [
                'specimen', 'assignment', 'assignment_id', 'query', 'subject', 'pident', 'coverage', 'ambig_count',
                'target_rank_id', 'target_rank_name', 'accession', 'tax_id', 'tax_name', 'target_rank', 'hi', 'low'
                ]),
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
    parser.add_argument('--max-ambiguous',
            default = 3,
            type = int,
            help = 'Maximum ambiguous count in reference sequences [%(default)s]')
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
    parser.add_argument('--group-label',
            default = 'all',
            help = 'Single group label for reads')
    parser.add_argument('--min-cluster-size',
            default = 0,
            type = int,
            help = 'minimum cluster size to include in classification output')
    parser.add_argument('--target-rank',
            help = 'Rank at which to classify. Default: "%(default)s"',
            default = 'species')
    parser.add_argument('--details-identity',
            help = 'Minimum identity threshold for classification details file',
            type = float,
            default = 90)
    parser.add_argument('--copy-numbers',
            type = Csv2Dict('tax_id', 'median'),
            default = {},
            help = '16S copy-number csv for correcting read numbers')
    parser.add_argument('--target-max-group-size',
            default = 3,
            type = int,
            help = """group multiple target-rank assignments that
                      excede a threshold to a higher rank [%(default)s]""")
    parser.add_argument('--group-def',
            action = 'append',
            default = [],
            help = """define a group threshold for a particular rank overriding --target-max-group-size.
                      example: genus:2""")

def get_copy_counts(taxids, copy_numbers, taxonomy, ranks):
    copy_counts = {}

    for t in (taxonomy[i] for i in taxids):
        tax_id = t['tax_id']

        if tax_id in copy_numbers:
            copy_counts[tax_id] = float(copy_numbers[tax_id])
        else:
            # return the copy_number for the lowest rank tax id available
            for r in ranks:
               if t[r] in copy_numbers:
                    copy_counts[tax_id] = float(copy_numbers[t[r]])
                    break

    return copy_counts

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
    condensed = (condense(g, floor_rank, max_size, ranks, rank_thresholds, target_rank) for _,g in groups)

    # flatten condensed
    condensed = [c for g in condensed for c in g]

    return condensed

def up_rank(tax, rank, ranks):
    """
    return the lowest possible defined rank id
    """

    index = ranks.index(rank)

    for r in ranks[:index]:
        if tax[r]:
            return tax[r]

    return '1' # root

def action(args):
    ### Rows
    etc = 'no match' # This row will hold all unmatched

    # groups have list position prioritization
    groups = [
        ('> {}%'.format(args.max_identity), lambda h: h['pident'] > args.max_identity),
        (None, lambda h: args.max_identity >= h['pident'] > args.min_identity),
        ('<= {}%'.format(args.min_identity), lambda h: h['pident'] <= args.min_identity),
    ]
    group_cats = map(itemgetter(0), groups)

    # assignment rank thresholds
    rank_thresholds = (d.split(':') for d in args.group_def)
    rank_thresholds = dict((k, int(v)) for k,v in rank_thresholds)

    ranks_rev = list(reversed(sequtils.RANKS))

    ### filter and format format blast data
    blast_results = DictReader(args.blast_file, fieldnames = sequtils.BLAST_HEADER)

    # some raw filtering
    blast_results = (dict(b, pident = float(b['pident'])) for b in blast_results)
    blast_results = (b for b in blast_results if b['pident'] >= args.details_identity)
    blast_results = (b for b in blast_results
            if float(args.weights.get(b['qseqid'], 1)) >= args.min_cluster_size)

    blast_results = (dict(b, coverage = coverage(b['qstart'], b['qend'], b['qlen'])) for b in blast_results)
    blast_results = (b for b in blast_results if b['coverage'] >= args.coverage)

    # add seq info
    blast_results = (dict(args.seq_info[b['sseqid']], **b) for b in blast_results)

    # add tax info
    blast_results = (dict(args.taxonomy[b['tax_id']], **b) for b in blast_results)

    # custom exclusion of tax records
    blast_results = (b for b in blast_results if b[args.target_rank] not in args.exclude_by_taxid)

    # (Optional) In some cases you want to filter some target hits
    blast_results = (b for b in blast_results if int(b['ambig_count']) <= args.max_ambiguous)
#    blast_results = (b for b in blast_results if b['qseqid'] != b['sseqid'])
    ###

    # add initial target rank id information
    blast_results = (dict(b, target_rank_id = up_rank(b, args.target_rank, ranks_rev)) for b in blast_results)

    # first, group by specimen
    if args.map:
        specimen_grouper = lambda s: args.map[s['qseqid']]
    elif args.all_one_group:
        specimen_grouper = lambda s: args.group_label
    else:
        specimen_grouper = lambda s: s['qseqid']

    blast_results = groupbyl(blast_results, key = specimen_grouper)

    assignments = [] # assignment list for assignment ids

    ### Columns
    out = DictWriter(args.out, extrasaction = 'ignore', fieldnames = [
        'max_percent', 'min_percent', 'max_coverage', 'min_coverage', 'assignment_id',
        'assignment', 'specimen', 'clusters', 'reads', 'pct_reads',
        'corrected', 'pct_corrected', 'target_rank', 'hi', 'low', 'tax_ids'
        ])
    out.writeheader()

    if args.out_detail:
        args.out_detail.writeheader()

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
                    target_ids = map(itemgetter('target_rank_id'), queries)
                    target_ids = frozenset(target_ids)

                    target_cats[target_ids].extend(queries)

                categories = dict(categories, **target_cats)

            # add query ids that were matched to a filter
            clusters |= set(map(itemgetter('qseqid'), matches))

            # remove all hits corresponding to a matched query id (cluster)
            hits = filter(lambda h: h['qseqid'] not in clusters, hits)

        # remaining hits go in the etc ('no match') category
        categories[etc] = hits

        # FIXME: will need to do this after grouping
        # needs taxids for 16S copy number corrections
        taxids = set(h['target_rank_id'] for v in categories.values() for h in v)

        copy_counts = get_copy_counts(taxids,
                                      args.copy_numbers,
                                      args.taxonomy,
                                      ranks_rev)

        # calculate read counts
        read_counts = ((t, set(map(itemgetter('qseqid'), h))) for t,h in categories.items())
        read_counts = ((t, sum(float(args.weights.get(q, 1)) for q in qs)) for t,qs in read_counts)
        read_counts = dict(read_counts)

        # corrected counts based on read_counts / mean(copy_counts)
        corrected_counts = ((t, set(map(itemgetter('tax_id'), h))) for t,h in categories.items())
        corrected_counts = ((c, mean(copy_counts.get(t, 1) for t in ts)) for c,ts in corrected_counts)
        corrected_counts = ((c, read_counts[c] / m if m else 0) for c,m in corrected_counts)
        corrected_counts = dict(corrected_counts)

        total_reads = sum(v for v in read_counts.values())
        total_corrected = sum(v for v in corrected_counts.values())

        # Print classifications per specimen sorted by # of reads in reverse (descending) order
        sort_by_reads = lambda (c,_): read_counts[c]
        for target_ids, hits in sorted(categories.items(), key = sort_by_reads, reverse = True):
            # only output categories with hits
            if hits:
                if target_ids not in assignments:
                    assignments.append(target_ids)

                assignment_id = assignments.index(target_ids)

                # assignment indexing
                taxids = set(map(itemgetter('tax_id'), hits))
                clusters = set(map(itemgetter('qseqid'), hits))
                coverages = set(map(itemgetter('coverage'), hits))
                percents = set(map(itemgetter('pident'), hits))

                reads = read_counts[target_ids]
                reads_corrected = corrected_counts[target_ids]

                # build tax name
                if target_ids in group_cats:
                    assignment = target_ids
                else:
                    names = [args.taxonomy[h['target_rank_id']]['tax_name'] for h in  hits]
                    selectors = [h['pident'] >= args.asterisk for h in hits]
                    assignment = sequtils.format_taxonomy(names, selectors, '*')

                out.writerow(dict(
                    hi = args.max_identity,
                    low = args.min_identity,
                    target_rank = args.target_rank,
                    specimen = specimen,
                    assignment_id = assignment_id,
                    assignment = assignment,
                    reads = int(reads),
                    pct_reads = '{0:.2f}'.format(reads / total_reads * 100),
                    corrected = int(reads_corrected),
                    pct_corrected = '{0:.2f}'.format(reads_corrected / total_corrected * 100),
                    clusters = len(clusters),
                    max_percent = '{0:.2f}'.format(max(percents)),
                    min_percent = '{0:.2f}'.format(min(percents)),
                    max_coverage = '{0:.2f}'.format(max(coverages)),
                    min_coverage = '{0:.2f}'.format(min(coverages)),
                    tax_ids = ' '.join(taxids)
                    ))

                if args.out_detail:
                    for h in hits:
                        args.out_detail.writerow(dict(
                        specimen = specimen,
                        assignment = assignment,
                        assignment_id = assignment_id,
                        hi = args.max_identity,
                        low = args.min_identity,
                        target_rank = args.target_rank,
                        **h))

