"""Classify sequences by grouping blast output by matching taxonomic names

Optional grouping by specimen and query sequences

The output is a table containing the following columns:

specimen, max_percent, min_percent, max_coverage, min_coverage,
assignment_id, assignment, clusters, reads, pct_read, corrected,
pct_corrected, target_rank, hi, low, tax_ids
"""

import sys
import logging

from os import path

import pandas as pd
import math

from bioy_pkg import sequtils
from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument(
        'blast_file', help='CSV tabular blast file of query and subject hits.')
    parser.add_argument(
        'seq_info', metavar='CSV',
        help='File mapping reference seq name to tax_id')
    parser.add_argument(
        'taxonomy', metavar='CSV',
        help="""Table defining the taxonomy for each tax_id""")
    parser.add_argument(
        '--copy-numbers', metavar='CSV',
        help="""Estimated 16s rRNA gene copy number for each tax_ids
        (CSV file with columns: tax_id, median)""")
    parser.add_argument(
        '--details-full', action='store_true',
        help='do not limit out_details to only larget cluster per assignment')
    parser.add_argument(
        '--group-def', metavar='INTEGER', action='append',
        default=[], help="""define a group threshold for a
        particular rank overriding --target-max-group-size. example:
        genus:2 (NOT IMPLEMENTED)""")
    parser.add_argument(
        '--has-header', action='store_true',
        help='specify this if blast data has a header')
    parser.add_argument(
        '--min-identity', default=90.0, metavar='PERCENT', type=float,
        help="""minimum identity threshold
        for accepting matches [>= %(default)s]""")
    parser.add_argument(
        '--max-identity', default=100.0, metavar='PERCENT', type=float,
        help="""maximum identity threshold for
        accepting matches [<= %(default)s]""")
    parser.add_argument(
        '--min-cluster-size', default=1, metavar='INTEGER', type=int,
        help="""minimum cluster size to include in
        classification output [%(default)s]""")
    parser.add_argument(
        '--min-coverage', default=95.0, type=float, metavar='PERCENT',
        help='percent of alignment coverage of blast result [%(default)s]')
    parser.add_argument(
        '-o', '--out', default=sys.stdout, type=Opener('w'),
        metavar='FILE',
        help="Classification results.")
    parser.add_argument(
        '-O', '--details-out', type=Opener('w'), metavar='FILE',
        help="""Optional details of taxonomic assignments.""")
    parser.add_argument(
        '--rank-thresholds', metavar='CSV',
        help="""Overrides target-rank.  Columns [tax_id,hi,low,rank]""")
    parser.add_argument(
        '--specimen', metavar='LABEL',
        help="""Single group label for reads""")
    parser.add_argument(
        '--specimen-map', metavar='CSV',
        help="""CSV file with columns (name, specimen) assigning sequences to
        groups. The default behavior is to treat all query sequences
        as belonging to one specimen.""")
    parser.add_argument(
        '--starred', default=100.0, metavar='PERCENT', type=float,
        help="""Names of organisms for which at least one reference
        sequence has pairwise identity with a query sequence of at
        least PERCENT will be marked with an asterisk [%(default)s]""")
    parser.add_argument(
        '--target-max-group-size', metavar='INTEGER', default=3, type=int,
        help="""group multiple target-rank assignments that excede a
        threshold to a higher rank [%(default)s]""")
    parser.add_argument(
        '--target-rank', default='species',
        help='Rank at which to classify. Default: "%(default)s"')
    parser.add_argument(
        '-w', '--weights', metavar='CSV',
        help="""Optional headless csv file with columns 'seqname',
        'count' providing weights for each query sequence described in
        the blast input (used, for example, to describe cluster sizes
        for corresponding cluster centroids).""")


def read_csv(filename, compression=None, **kwargs):
    """read a csv file using pandas.read_csv with compression defined by
    the file suffix unless provided.

    """

    suffixes = {'.bz2': 'bz2', '.gz': 'gzip'}
    compression = compression or suffixes.get(path.splitext(filename)[-1])
    kwargs['compression'] = compression

    return pd.read_csv(filename, **kwargs)


def action(args):
    # pd.set_option('display.max_columns', None)
    # pd.set_option('display.max_rows', None)

    # format blast data and add additional available information
    names = None if args.has_header else sequtils.BLAST_HEADER
    header = 0 if args.has_header else None
    usecols = ['qseqid', 'sseqid', 'pident', 'coverage']
    blast_results = read_csv(
        args.blast_file,
        names=names,
        na_filter=True,  # False is faster
        header=header,
        usecols=usecols)

    # get a set of qseqids for identifying [no blast hits] after filtering
    qseqids = blast_results[['qseqid']].drop_duplicates().set_index('qseqid')

    # filter out low coverage hits
    blast_results = blast_results[
        blast_results['coverage'] >= args.min_coverage]

    seq_info = read_csv(
        args.seq_info,
        usecols=['seqname', 'tax_id', 'accession'],
        dtype=dict(tax_id=str),
        index_col='seqname')
    # rename index to match blast results
    seq_info.index.name = 'sseqid'

    # merge blast results with seq_info - do this early so that
    # refseqs not represented in the blast results are discarded in
    # the merge.
    blast_results = blast_results.join(seq_info, on='sseqid')

    # load taxonomy
    tax_cols = ['tax_id', 'tax_name']
    if args.rank_thresholds:
        rank_thresholds = read_csv(
            args.rank_thresholds,
            dtype=dict(tax_id=str))
        rank_thresholds.index = rank_thresholds['tax_rank'].apply(
            lambda x: sequtils.RANKS.index(x))
        # sort items by RANK specificity
        rank_thresholds = rank_thresholds.sort(ascending=False)
        rank_thresholds = rank_thresholds.reset_index(drop=True)

        # we will need tax columns from both the target_rank and tax_rank
        ranks = set(rank_thresholds['target_rank'].unique())

        # calculate the floor rank from most specific target_rank
        floor_rank = sorted(
            ranks, key=lambda x: sequtils.RANKS.index(x), reverse=True)[0]

        # get the tax_ranks as well and combine them with the target_ranks
        tax_ranks = set(rank_thresholds['tax_rank'].unique())
        ranks = list(ranks | tax_ranks)
    else:
        floor_rank = args.target_rank
        ranks = [args.target_rank]

    # TODO: consider getting rank orderding from columns of taxomomy table
    taxonomy = read_csv(args.taxonomy, dtype=str)
    # set index after assigning dtype and preserve tax_id column for later
    taxonomy = taxonomy.set_index('tax_id', drop=False)  # keep tax_id column

    blast_results = blast_results.join(taxonomy[ranks], on='tax_id')
    #  tax_id will be later set from the target_rank
    blast_results = blast_results.drop('tax_id', axis=1)


# - construct table hits by joining blast_results and seq_info (seqname, tax_id, pident); add boolean column classified
# - construct table tax_info by joining seq_info with taxonomy (seqname, tax_id, column for each rank) [may need to fill in missing values for some or all ranks]

# for rank, threshold in thresholds:
# - create a table by:
#   1. join rows of hits where classified is False and pident > threshold (could use a second table to look up tax-id specific thresholds at this rank) with tax_info using tax_id --> seqname, pident, target_tank, tax_id (tax_id at this target rank corresponding to the appropriate column in tax_info)
#   2. set hits.classified = True for all seqnames represented in the table

# stack up these tables, group by seqname, and combine names

    # assign target rank
    if args.rank_thresholds:
        target_ranks = []
        # iter through sorted/prioritized rank_thresholds
        for i, threshold in rank_thresholds.iterrows():
            # filter for hits that match the threshold
            # TODO: is it ok to use <= (vs <) in the lower range?
            matches = blast_results[
                (threshold['tax_id'] == blast_results[threshold['tax_rank']]) &
                (threshold['hi'] >= blast_results['pident']) &
                (threshold['low'] <= blast_results['pident'])]

            # pivot tax_rank and tax_id and drop for merging
            threshold[threshold['tax_rank']] = threshold['tax_id']
            threshold = threshold.drop(['tax_rank', 'tax_id'])
            threshold = pd.DataFrame(threshold.to_dict(), index=[i])
            threshold['target_priority'] = i

            # append merged matches to target ranks and remove from
            # blast_results
            target_ranks.append(matches.merge(threshold))
            not_matches = blast_results.index.diff(matches.index)
            blast_results = blast_results.loc[not_matches]

        # put blast hits back together with specified target ranks
        blast_results = pd.concat(target_ranks)

        # FIXME:CR - this is very time consuming, is there a better way?
        def slim_results(df):
            """
            Remove less specific target_rank hits.

            If a hit's target rank shares taxonomy with a more specific
            target hit then it will be dropped in favor
            of the more specific target hit(s)
            """
            slimmed = pd.DataFrame(columns=df.columns.tolist() + tax_cols)
            slimmed['target_priority'] = slimmed['target_priority'].astype(int)
            df = df.groupby(by=['target_priority', 'target_rank'])
            for (_, target_rank), vals in df:
                curated_ids = slimmed[target_rank].drop_duplicates()
                vals = vals[~vals[target_rank].isin(curated_ids)]
                vals = vals.join(taxonomy[tax_cols], on=target_rank)
                slimmed = slimmed.append(vals)
            return slimmed

        blast_results = blast_results.groupby(
            by='qseqid', sort=False, group_keys=False).apply(slim_results)
    else:
        blast_results = blast_results[
            blast_results['pident'] <= args.max_identity]
        blast_results = blast_results[
            blast_results['pident'] >= args.min_identity]
        blast_results = blast_results.join(
            taxonomy[tax_cols], on=args.target_rank)
        blast_results['target_rank'] = args.target_rank
        blast_results['target_priority'] = 0
        blast_results['hi'] = args.max_identity
        blast_results['low'] = args.min_identity

    # TODO: up-rank hits that happen to have no tax_id at the given target rank
    blast_results = blast_results[blast_results['tax_id'].notnull()]

    # merge qseqids that have no hits back into blast_results
    blast_results = qseqids.join(
        blast_results.set_index('qseqid'), how='outer').reset_index()

    # assign specimen groups
    specimens = blast_results[['qseqid']].drop_duplicates()
    specimens = specimens.set_index('qseqid')

    # load specimen-map and assign specimen names
    if args.specimen_map:
        spec_map = read_csv(
            args.specimen_map,
            names=['qseqid', 'specimen'],
            index_col='qseqid')
        specimens = specimens.join(spec_map)
    elif args.specimen:
        specimens['specimen'] = args.specimen
    else:
        specimens['specimen'] = specimens.index

    # join back into results.  Now all hits have an assigned specimen label
    blast_results = blast_results.join(specimens, on='qseqid')

    # assign seqs that had no results to [no blast_result]
    no_hits = blast_results[blast_results.sseqid.isnull()]
    no_hits['assignment'] = '[no blast result]'
    no_hits['assignment_hash'] = 0

    # move on to seqs that have blast hits
    blast_results = blast_results[blast_results.sseqid.notnull()]

    # TODO: this is slow, integrate pandas into sequtils.condense_ids
    tax_dict = {i: t.to_dict() for i, t in taxonomy.fillna('').iterrows()}

    # create mapping from tax_id to its condensed id and set assignment hash
    def condense_ids(df):
        condensed = sequtils.condense_ids(
            df.tax_id.unique(),
            tax_dict,
            floor_rank=floor_rank,
            max_size=args.target_max_group_size)
        condensed = pd.DataFrame(
            condensed.items(),
            columns=['tax_id', 'condensed_id'])
        condensed = condensed.set_index('tax_id')
        assignment_hash = hash(frozenset(condensed.condensed_id.unique()))
        condensed['assignment_hash'] = assignment_hash
        return df.join(condensed, on='tax_id')

    # create condensed assignment hashes by qseqid
    blast_results = blast_results.groupby(
        by=['qseqid'], sort=False).apply(condense_ids)

    # star hits if one hit meets star threshold
    def star(df):
        df['starred'] = df.pident.apply(lambda x: x >= args.starred).any()
        return df

    # star condensed ids
    blast_results = blast_results.groupby(
        by=['assignment_hash', 'condensed_id'], sort=False).apply(star)

    def assign(df):
        ids_stars = df.groupby(by=['condensed_id', 'starred']).groups.keys()
        df['assignment'] = sequtils.compound_assignment(ids_stars, tax_dict)
        return df

    blast_results = blast_results.groupby(
        by=['assignment_hash'], sort=False).apply(assign)

    # put assignments and no assignments back together
    blast_results = pd.concat([blast_results, no_hits])

    # concludes our blast details, on to summarizing output

    # now for some assignment grouping and summarizing

    def agg_columns(df):
        agg = {}
        agg['hi'] = df['hi'].max()
        agg['low'] = df['low'].min()
        agg['target_rank'] = df.sort(
            columns='target_priority').head(1)['target_rank']
        agg['max_percent'] = df['pident'].max()
        agg['min_percent'] = df['pident'].min()
        return pd.DataFrame(agg)

    output = blast_results.set_index(
        ['specimen', 'assignment_hash', 'assignment'])
    output = output.groupby(
        level=['specimen', 'assignment_hash'],
        group_keys=False).apply(agg_columns)
    output = output.reset_index(level='assignment')

    # qseqid cluster stats
    clusters = blast_results[['qseqid', 'specimen', 'assignment_hash']]
    clusters = clusters.drop_duplicates().set_index('qseqid')

    if args.weights:
        weights = read_csv(
            args.weights,
            names=['qseqid', 'weight'],
            index_col='qseqid')
        clusters = clusters.join(weights)
        # switch back to int and set no info to weight of 1
        clusters['weight'] = clusters['weight'].fillna(1).astype(int)
    else:
        clusters['weight'] = 1

    clusters = clusters.groupby(by=['specimen', 'assignment_hash'], sort=False)

    output['reads'] = clusters['weight'].sum()

    def freq(df):
        df['pct_reads'] = df['reads'] / df['reads'].sum() * 100
        return df

    output = output.groupby(level='specimen').apply(freq)

    output['clusters'] = clusters.size()

    # copy number corrections
    if args.copy_numbers:
        copy_numbers = read_csv(
            args.copy_numbers,
            dtype=dict(tax_id=str),
            usecols=['tax_id', 'median']).set_index('tax_id')

        # get root out '1' and set it as the
        # default correction value with index nan
        default = copy_numbers.get_value('1', 'median')
        default = pd.DataFrame(default, index=[None], columns=['median'])
        copy_numbers = copy_numbers.append(default)

        # do our copy number correction math
        corrections = blast_results[['tax_id', 'specimen', 'assignment_hash']]
        corrections = corrections.drop_duplicates().set_index('tax_id')
        corrections = corrections.join(copy_numbers)
        corrections = corrections.groupby(
            by=['specimen', 'assignment_hash'], sort=False)
        corrections = corrections['median'].mean()
        output['corrected'] = output['reads'] / corrections

        def pct_corrected(df):
            df['pct_corrected'] = df['corrected'] / df['corrected'].sum() * 100
            return df

        output = output.groupby(level='specimen').apply(pct_corrected)

        # reset corrected counts to int
        output['corrected'] = output['corrected'].apply(math.ceil)
        output['corrected'] = output['corrected'].astype(int)

    # sort by read count and specimen
    columns = ['corrected'] if 'corrected' in output else ['reads']
    output = output.sort(columns=columns, ascending=False)
    output = output.reset_index(level='assignment_hash', drop=True)
    output = output.sort_index()

    # create assignment ids by specimen
    def assignment_id(df):
        df = df.reset_index(drop=True)  # specimen is retained in the group key
        df.index.name = 'assignment_id'
        return df

    output = output.groupby(level="specimen", sort=False).apply(assignment_id)

    # output results
    with args.out as out:
        output.to_csv(out, index=True, float_format='%.2f')

    # output to details.csv.bz2
    if args.details_out:
        blast_results = blast_results.merge(output.reset_index(), how='left')

        if not args.details_full:
            largest = clusters.apply(lambda x: x['weight'].nlargest(1))
            blast_results = blast_results.merge(largest.reset_index())

        # nobody wants to see the assignment hash
        blast_results = blast_results.drop(labels='assignment_hash', axis=1)

        with args.details_out as out_details:
            blast_results.to_csv(
                out_details, header=True, index=False, float_format='%.2f')
