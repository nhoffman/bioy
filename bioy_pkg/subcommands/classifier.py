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

    # input files
    parser.add_argument(
        'blast_file', help='CSV tabular blast file of query and subject hits.')
    parser.add_argument(
        '-w', '--weights', metavar='CSV',
        help="""Optional headless csv file with columns 'seqname',
        'count' providing weights for each query sequence described in
        the blast input (used, for example, to describe cluster sizes
        for corresponding cluster centroids).""")
    parser.add_argument(
        '-s', '--seq-info', required=True, metavar='CSV',
        help='File mapping reference seq name to tax_id')
    parser.add_argument(
        '-t', '--taxonomy', required=True, metavar='CSV',
        help="""Table defining the taxonomy for each tax_id""")
    parser.add_argument(
        '--copy-numbers', metavar='CSV',
        help="""Estimated 16s rRNA gene copy number for each tax_ids
        (CSV file with columns: tax_id, median)""")
    parser.add_argument(
        '-m', '--specimen-map', metavar='CSV',
        help="""CSV file with columns (name, specimen) assigning sequences to
        groups. The default behavior is to treat all query sequences
        as belonging to one specimen.""")
    # input file parameters
    parser.add_argument(
        '--has-header', action='store_true', default=False,
        help='specify this if blast data has a header')
    # output files
    parser.add_argument(
        '-o', '--out', default=sys.stdout, type=Opener('w'),
        metavar='FILE',
        help="Classification results.")
    parser.add_argument(
        '-O', '--out-details', type=Opener('w'), metavar='FILE',
        help="""Optional details of taxonomic assignments.""")
    # classification parameters
    parser.add_argument(
        '--rank', default='species',
        help='Rank at which to classify. Default: "%(default)s"')
    parser.add_argument(
        '--min-identity', default=90.0, metavar='PERCENT', type=float,
        help="""minimum identity threshold
        for accepting matches [> %(default)s]""")
    parser.add_argument(
        '--max-identity', default=100.0, metavar='PERCENT', type=float,
        help="""maximum identity threshold for
        accepting matches [<= %(default)s]""")
    parser.add_argument(
        '--starred', default=100.0, metavar='PERCENT', type=float,
        help="""Names of organisms for which at least one reference
        sequence has pairwise identity with a query sequence of at
        least PERCENT will be marked with an asterisk [%(default)s]""")
    parser.add_argument(
        '--min_coverage', default=95.0, type=float, metavar='PERCENT',
        help='percent of alignment coverage of blast result [%(default)s]')
    parser.add_argument(
        '--min-cluster-size', default=1, metavar='INTEGER', type=int,
        help='minimum cluster size to include in classification output')
    parser.add_argument(
        '--all-one-group', action='store_true', default=False,
        help="""If --map is not provided, the default behavior is to
        treat all reads as one group; use this option to treat each
        read as a separate group [%(default)s]""")

    # options for behavior of naming aglorithm
    parser.add_argument(
        '--target-max-group-size', metavar='INTEGER', default=3, type=int,
        help="""group multiple target-rank assignments that excede a
        threshold to a higher rank [%(default)s]""")
    parser.add_argument(
        '--group-def', metavar='INTEGER', action='append',
        default=[], help="""define a group threshold for a
        particular rank overriding --target-max-group-size. example:
        genus:2""")

    # TODO: what does this do?
    # CR: text label to handle all blast results as one group when no specimen
    # map is specified.  Goes with the all-one-group option above
    parser.add_argument(
        '--group-label', metavar='LABEL', default='all',
        help="""Single group label for reads""")

    parser.add_argument(
        '--details-full', action='store_true',
        help='do not limit out_details to only larget cluster per assignment')


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

    seq_info = read_csv(
        args.seq_info,
        usecols=['seqname', 'tax_id', 'accession'],
        dtype=dict(tax_id=str,
                   length=int,
                   ambig_count=int,
                   is_type=bool),
        index_col='seqname')
    seq_info.index.name = 'sseqid'

    # must set index after assigning datatype
    taxonomy = read_csv(args.taxonomy, dtype=str).set_index('tax_id')

    # format blast data and add additional available information
    names = None if args.has_header else sequtils.BLAST_HEADER[:6]
    header = 0 if args.has_header else None
    blast_results = read_csv(args.blast_file,
                             names=names,
                             na_filter=True,  # False (faster)
                             header=header,
                             usecols=['qseqid', 'sseqid', 'pident'])

    # filter out high and low identity hits and low coverage hits
    blast_results = blast_results[blast_results['pident'] >= args.min_identity]
    blast_results = blast_results[blast_results['pident'] <= args.max_identity]

    if 'coverage' in blast_results:
        blast_results = blast_results[
            blast_results['coverage'] >= args.min_coverage]

    # hi and low
    blast_results['hi'] = args.max_identity
    blast_results['low'] = args.min_identity

    # assign target rank
    # TODO: read in thresholds file
    blast_results['target_rank'] = args.rank

    # merge blast results with seq_info - do this early so that
    # refseqs not represented in the blast results are discarded in
    # the merge.
    blast_results = blast_results.join(seq_info, on='sseqid')

    # merge with taxonomy to associate original tax_ids with tax_ids
    # at the specified rank.
    blast_results = blast_results.join(taxonomy[[args.rank]], on='tax_id')

    # merge with taxonomy again to get tax_names
    blast_results = blast_results.join(taxonomy[['tax_name']], on=args.rank)

    # assign seqs that had no results to [no blast_result]
    no_hits = blast_results[blast_results.sseqid.isnull()]
    no_hits['assignment'] = '[no blast result]'
    no_hits['assignment_hash'] = 0

    # move on to seqs that had hits
    blast_results = blast_results[blast_results.sseqid.notnull()]

    # TODO: this is taking a lot of time, integrate pandas into sequtils
    tax_dict = {i: t.to_dict() for i, t in taxonomy.fillna('').iterrows()}

    # create mapping from tax_id to its condensed id and set assignment hash
    def condense_ids(df):
        condensed = sequtils.condense_ids(df.tax_id.unique(), tax_dict)
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

    # concludes our blast details

    # now for some assignment grouping and summarizing
    output = blast_results[['assignment_hash', 'assignment', 'target_rank']]
    output = output.drop_duplicates().set_index('assignment_hash')

    # the qseqid cluster stats
    qseqids = blast_results[['qseqid', 'assignment_hash']]
    qseqids = qseqids.drop_duplicates().set_index('qseqid')

    if args.weights:
        weights = read_csv(args.weights,
                           header=None,
                           names=['qseqid', 'weight'],
                           index_col='qseqid')
        qseqids = qseqids.join(weights)
        # switch back to int and set no info to weight of 1
        qseqids['weight'] = qseqids['weight'].fillna(1).astype(int)
    else:
        qseqids['weight'] = 1

    qseqids = qseqids.groupby(by=['assignment_hash'], sort=False)

    output['reads'] = qseqids.weight.sum()
    output['freq'] = output['reads'] / output['reads'].sum() * 100
    output['clusters'] = qseqids.size()

    if args.copy_numbers:
        copy_numbers = read_csv(
            args.copy_numbers,
            dtype=dict(tax_id=str),
            usecols=['tax_id', 'median']).set_index('tax_id')

        # get root out '1' and set it as the default correction value
        default = copy_numbers.get_value('1', 'median')
        default = pd.DataFrame(default, index=[None], columns=['median'])
        copy_numbers = copy_numbers.append(default)

        # do our copy number correction math
        corrections = blast_results[['tax_id', 'assignment_hash']]
        corrections = corrections.drop_duplicates().set_index('tax_id')
        corrections = corrections.join(copy_numbers)
        corrections = corrections.groupby(by=['assignment_hash'], sort=False)
        corrections = corrections['median'].mean()
        output['corrected'] = output['reads'] / corrections
        pct_corrected = output['corrected'] / output['corrected'].sum() * 100
        output['pct_corrected'] = pct_corrected
        # reset corrected counts to int
        output['corrected'] = output['corrected'].apply(math.ceil)
        output['corrected'] = output['corrected'].astype(int)

    # qseqid hit ranges
    assignment_groups = blast_results.groupby(
        by=['assignment_hash'], sort=False)
    output['max_percent'] = assignment_groups['pident'].max()
    output['min_percent'] = assignment_groups['pident'].min()

    output = output.sort('corrected', ascending=False)
    output.to_csv(args.out, index=False, float_format='%.2f')

    # output to details.csv.bz2
    if args.out_details:
        with args.out_details as out_details:
            if args.details_full:
                blast_results.to_csv(out_details,
                                     header=True,
                                     index=False,
                                     float_format='%.2f')
            else:
                largest = qseqids.apply(lambda x: x['weight'].nlargest(1))
                largest = largest.reset_index()
                largest.merge(blast_results).to_csv(out_details,
                                                    header=True,
                                                    index=False,
                                                    float_format='%.2f')
