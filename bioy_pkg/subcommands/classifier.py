"""Classify sequences by grouping blast output by matching taxonomic names

Optional grouping by specimen and query sequences

TODO: describe the algorithm here...

 Input details
===============

rank_thresholds
---------------

This table defines the similarity thresholds at rank. The structure is
as follows::

  +--------+-----+--------------|
  | tax_id | low | target_rank  |
  +--------+-----+--------------|
  | 1      | 99  | species      |
  | 1      | 97  | genus        |
  | 1      | 95  | family       |
  +--------+-----+--------------|

The `tax_id` column identifies the subtree of the taxonomy to which
the threshold defined in `low` should be applied. In the example
above, the threshold of 99 applies to all tax_ids with an ancestor of
tax_id=1 (ie, the entire taxonomy). A more useful example might be this::

 |--------+-----+-------------|
 | tax_id | low | target_rank |
 |--------+-----+-------------|
 |   2049 |  97 | species     |
 |   2049 |  95 | genus       |
 |   2049 |  93 | family      |
 |--------+-----+-------------|

Here, all organisms belonging to family Actinomycetaceae (tax_id 2049)
will be classified to the species level using a threshold of 97%
instead of 99% in the first example, because we know that there is
more species-level heterogeneity within this family.

seq_info
--------

A csv file with columns "seqname","tax_id" and an optional column
"accession", which gets passed through into the details output.


Output
======

* `--out` is a table containing the following columns:

specimen, max_percent, min_percent, max_coverage, min_coverage,
assignment_id, assignment, clusters, reads, pct_read, corrected,
pct_corrected, target_rank, low, tax_ids



"""

import sys
import logging

from os import path

import pandas as pd
import math

from bioy_pkg import sequtils
from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)

def freq(df):
    df['pct_reads'] = df['reads'] / df['reads'].sum() * 100
    return df

def read_csv(filename, compression=None, **kwargs):
    """read a csv file using pandas.read_csv with compression defined by
    the file suffix unless provided.

    """

    suffixes = {'.bz2': 'bz2', '.gz': 'gzip'}
    compression = compression or suffixes.get(path.splitext(filename)[-1])
    kwargs['compression'] = compression

    return pd.read_csv(filename, **kwargs)

def slim_results(df):
    """Drop unqualified hits and return remaining hits above the identiy
    threshold defined for this rank.

    """

    # first drop all rows corresponding to hits that don't meet the
    # threshold at the specified target rank
    df = df[df['low'] <= df['pident']]

    # among the remaining hits, keep hits at the most stringent
    # threshold.
    return df[df['low'].max() <= df['pident']]


def star(df, starred):
    df['starred'] = df.pident.apply(lambda x: x >= starred).any()
    return df


def condense_ids(df, tax_dict, max_group_size):

    """ create mapping from tax_id to its
        condensed id and set assignment hash """
    condensed = sequtils.condense_ids(
        df.tax_id.unique(),
        tax_dict,
        max_size=max_group_size)
    condensed = pd.DataFrame(
        condensed.items(),
        columns=['tax_id', 'condensed_id']).set_index('tax_id')
    assignment_hash = hash(frozenset(condensed.condensed_id.unique()))
    condensed['assignment_hash'] = assignment_hash
    return df.join(condensed, on='tax_id')


def assign(df, tax_dict):
    ids_stars = df.groupby(by=['condensed_id', 'starred']).groups.keys()
    df['assignment'] = sequtils.compound_assignment(ids_stars, tax_dict)
    return df


def agg_columns(df):
    agg = {}
    agg['low'] = df['low'].min()
    agg['max_percent'] = df['pident'].max()
    agg['min_percent'] = df['pident'].min()
    target_rank = df['target_rank'].dropna().drop_duplicates()
    if target_rank.empty:
        target_rank = None
    else:
        specificity = lambda x: sequtils.RANKS.index(x)
        target_rank.index = target_rank.apply(specificity)
        target_rank = target_rank.sort_index()
        target_rank = target_rank.iloc[-1]
    agg['target_rank'] = target_rank
    return pd.Series(agg)


def pct_corrected(df):
    df['pct_corrected'] = df['corrected'] / df['corrected'].sum() * 100
    return df

def assignment_id(df):
    df = df.reset_index(drop=True)  # specimen is retained in the group key
    df.index.name = 'assignment_id'
    return df

# TODO: organize options as 1) all inputs 2) all outputs 3) other options
def build_parser(parser):
    parser.add_argument(
        'blast_file', metavar='blast.csv(.bz2)',
        help='CSV tabular blast file of query and subject hits.')
    parser.add_argument(
        'seq_info', metavar='seq_info.csv(.bz2)',
        help='File mapping reference seq name to tax_id')
    parser.add_argument(
        'taxonomy', metavar='taxonomy.csv(.bz2)',
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
        '--min-identity', default=0.0, metavar='PERCENT', type=float,
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
        '--min-coverage', default=0.0, type=float, metavar='PERCENT',
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
        help="""Columns [tax_rank,tax_id,low,rank]""")
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


def action(args):
    # format blast data and add additional available information
    names = None if args.has_header else sequtils.BLAST_HEADER
    header = 0 if args.has_header else None
    usecols = ['qseqid', 'sseqid', 'pident', 'coverage']
    blast_results = read_csv(
        args.blast_file,
        dtype=dict(qseqid=str, sseqid=str, pident=float, coverage=float),
        names=names,
        na_filter=True,  # False is faster
        header=header,
        usecols=usecols)

    # get a set of qseqids for identifying [no blast hits] after filtering
    qseqids = blast_results[['qseqid']].drop_duplicates().set_index('qseqid')

    # filter out low coverage hits
    if args.min_coverage:
        blast_results = blast_results[
            blast_results['coverage'] >= args.min_coverage]

    # Define default thresholds; these will be overwritten by rank-specific thresholds below.
    # apply pident ceiling
    blast_results = blast_results[
        blast_results['pident'] <= args.max_identity]

    # apply floor
    blast_results = blast_results[
        blast_results['pident'] >= args.min_identity]

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

    # TODO: consider getting rank orderding from columns of taxomomy table
    usecols = ['tax_id', 'tax_name', 'rank'] + sequtils.RANKS
    taxonomy = read_csv(args.taxonomy, usecols=usecols, dtype=str)
    # set index after assigning dtype and preserve tax_id column
    taxonomy = taxonomy.set_index('tax_id')

    targeted = []
    if args.rank_thresholds:
        # Each tax_id T in the tax_id column of rank_thresholds
        # identifies a subtree of the taxonomy. The value in "low"
        # provides the threshold for tax_id T; all child tax_ids in
        # this subtree should be given a value of NULL.

        rank_thresholds = read_csv(
            args.rank_thresholds,
            dtype=dict(tax_id=str))

        # create a data_frame with columns tax_id, low defining the
        # classification threshold for each tax_id using
        # rank_thresholds.
        thresholds =
        # tax_id, low
        # '543', 95
        # '570', NA
        # '571', NA

        blast_results = blast_results.join(
            taxonomy[sequtils.RANKS], on='tax_id')

        # sort thresholds by rank and low threshold
        # rank_thresholds['rank_index'] = rank_thresholds['tax_rank'].apply(
        #     lambda x: sequtils.RANKS.index(x))
        # rank_thresholds = rank_thresholds.sort(
        #     columns=['rank_index', 'low'], ascending=False)

        # iterate over target_ranks, most specific first
        hitlist = []
        for target_rank in [rank for rank in sequtils.RANKS if rank in rank_thresholds.target_rank]:
            # define the threshold "low" for each tax_id at this target rank

            # use the tax_ids for this target_rank to define the
            # classification threshold for this iteration by joining
            # thresholds on thresholds.tax_id =
            # blast_results.<target_rank>
            # lows is a series with length = the number of rows in blast_results
            lows = <threshold defined in thresholds for the column in blast_results corresponding to the current target_rank>
            hits = blast_results[blast_results['pident'] >= lows]
            hits['target_rank'] = target_rank
            hitlist.append(hits)

            blast_results = (all rows in blast_results not in hits, but not including any qseqs represented in hits)

        # any qsesqs remaining in blast_results should be classified as [no blast result]

        # concatenate hitlist, group by qseq, and create assignments based on
        # remaining reference sequences for each qseq

        for i, threshold in rank_thresholds.iterrows():
            # threshold['tax_rank'] is the column containing the
            # tax_id (provided in threshold['tax_id']) to whose
            # children (in column 'target_rank') we will assign the
            # threshold ('low') defined in this iteration.

            hits = blast_results[
                (blast_results[threshold['tax_rank']] == threshold['tax_id']) &
                (blast_results['pident'] >= threshold['low'])]

            hits['low'] = threshold['low']
            hits['target_rank'] = threshold['target_rank']

            # retain tax_id of target rank only
            hits = hits.drop('tax_id', axis=1)

            # define tax_ids at target_rank
            hits['tax_id'] = hits[threshold['target_rank']]

            # keep hits not meeting the threshold for the next iteration
            diff = blast_results.index.diff(hits.index)
            blast_results = blast_results.loc[diff]
            targeted.append(hits)


    blast_results['low'] = args.min_identity
    blast_results['target_rank'] = args.target_rank
    blast_results = pd.concat(targeted + [blast_results])

    blast_results = blast_results.join(
        taxonomy[['tax_name', 'rank']], on='tax_id')

    # keep most specific hits
    blast_results = blast_results.groupby(
        by=['qseqid'], sort=False, group_keys=False).apply(slim_results)

    # merge qseqids that have no hits back into blast_results
    blast_results = blast_results.join(qseqids, on='qseqid', how='outer')

    # assign specimen groups
    specimens = blast_results[['qseqid']].drop_duplicates().set_index('qseqid')

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

    # create condensed assignment hashes by qseqid
    blast_results = blast_results.groupby(
        by=['qseqid'], sort=False).apply(lambda x: condense_ids(x, tax_dict, args.target_max_group_size))

    # star condensed ids if one hit meets star threshold
    blast_results = blast_results.groupby(
        by=['assignment_hash', 'condensed_id'], sort=False).apply(lambda x: star(x, args.starred))

    blast_results = blast_results.groupby(
        by=['assignment_hash'], sort=False).apply(lambda x: assign(x, tax_dict))

    # put assignments and no assignments back together
    blast_results = pd.concat([blast_results, no_hits])

    # concludes our blast details, on to assignments and output summary

    output = blast_results.groupby(
        by=['specimen', 'assignment_hash', 'assignment'])
    output = output.apply(agg_columns)
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
        # enforce weight dtype as int and unlisted qseq's to weight of 1
        clusters['weight'] = clusters['weight'].fillna(1).astype(int)
    else:
        clusters['weight'] = 1

    clusters = clusters.groupby(by=['specimen', 'assignment_hash'], sort=False)

    output['reads'] = clusters['weight'].sum()

    output = output.groupby(level='specimen').apply(freq)

    output['clusters'] = clusters.size()

    # copy number corrections
    if args.copy_numbers:
        copy_numbers = read_csv(
            args.copy_numbers,
            dtype=dict(tax_id=str),
            usecols=['tax_id', 'median']).set_index('tax_id')

        # get root out (taxid: 1) and set it as the
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

        # reset corrected counts to int before calculating pct_corrected
        output['corrected'] = output['corrected'].apply(math.ceil)
        output['corrected'] = output['corrected'].fillna(1).astype(int)

        output = output.groupby(level='specimen').apply(pct_corrected)

    # round up any pct < 0.01 for before sorting
    round_up = lambda x: max(0.01, x)
    output['pct_reads'] = output['pct_reads'].map(round_up)
    if args.copy_numbers:
        output['pct_corrected'] = output['pct_corrected'].map(round_up)

    # sort by:
    # 1) specimen
    # 2) read/corrected count
    # 3) cluster count
    # 4) alpha assignment
    columns = ['corrected'] if args.copy_numbers else ['reads']
    columns += ['clusters', 'assignment']
    output = output.sort(columns=columns, ascending=False)
    output = output.reset_index(level='assignment_hash', drop=True)
    output = output.sort_index()

    # create assignment ids by specimen
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

        # nobody cares about the assignment hash
        blast_results = blast_results.drop(labels='assignment_hash', axis=1)

        with args.details_out as out_details:
            blast_results.to_csv(
                out_details, header=True, index=False, float_format='%.2f')
