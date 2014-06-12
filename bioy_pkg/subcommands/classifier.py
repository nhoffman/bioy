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

import numpy as np
import pandas as pd

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
        '-g', '--groups', metavar='CSV',
        help="""CSV file with columns (name, specimen) assigning sequences to
        groups. The default behavior is to treat all query sequences
        as belonging to one specimen.""")
    parser.add_argument(
        '--exclude-taxids', metavar='File',
        help="""A file containing taxids - blast hits to reference sequences with
        these taxids will be discarded.""")

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
        '--min-identity', default=99.0, metavar='PERCENT', type=float,
        help='minimum identity threshold for accepting matches [> %(default)s]')
    parser.add_argument(
        '--max-identity', default=100.0, metavar='PERCENT', type=float,
        help='maximum identity threshold for accepting matches [<= %(default)s]')
    parser.add_argument(
        '--starred', default=100.0, metavar='PERCENT', type=float,
        help="""Names of organisms for which at least one reference
        sequence has pairwise identity with a query sequence of at
        least PERCENT will be marked with an asterisk [%(default)s]""")
    parser.add_argument(
        '-c', '--coverage', default=95.0, type=float, metavar='PERCENT',
        help='percent of alignment coverage of blast result [%(default)s]')
    parser.add_argument(
        '--min-cluster-size', default=1, metavar='INTEGER', type=int,
        help='minimum cluster size to include in classification output')
    parser.add_argument(
        '--max-ambiguous', metavar='INTEGER', default=3, type=int,
        help='Maximum ambiguous count in reference sequences [%(default)s]')
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

    # options for details file
    parser.add_argument(
        '--details-identity', metavar='PERCENT', type=float, default=90.0,
        help='Minimum identity to include blast hits in details file')
    parser.add_argument(
        '--details-full', action='store_true', default='False',
        help='do not limit out_details to only larget cluster per assignment')

def read_csv(filename, compression=None, **kwargs):
    """read a csv file using pandas.read_csv with compression defined by
    the file suffix unless provided.

    """

    suffixes = {'.bz2': 'bz2', '.gz': 'gzip'}
    compression = compression or suffixes.get(path.splitext(filename)[-1])
    kwargs['compression'] = compression

    return pd.read_csv(filename, **kwargs)


# details columns (maybe)
# ['specimen', 'assignment', 'assignment_id', 'qseqid', 'sseqid',
# 'pident', 'coverage', 'ambig_count', 'accession', 'tax_id',
# 'tax_name', 'target_rank', 'rank', 'hi', 'low']

def action(args):
#    pd.set_option('display.max_columns', None)
#    pd.set_option('display.max_rows', None)

    seq_info = read_csv(args.seq_info,
                        dtype={'tax_id':str,
                               'length': int,
                               'ambig_count': int,
                               'is_type': bool})
    seq_info = seq_info.rename(columns=dict(seqname='sseqid'))
    seq_info = seq_info.set_index('sseqid')

    taxonomy = read_csv(args.taxonomy, dtype=str)
    taxonomy = taxonomy.set_index('tax_id')

    # format blast data and add additional available information
    blast_results = read_csv(args.blast_file,
                             header=0 if args.has_header else None)

    if not args.has_header:
        columns = dict(zip(blast_results.columns, sequtils.BLAST_HEADER))
        blast_results = blast_results.rename(columns=columns)

    # add a column indicating if a hit meets the threshold for
    # starring.
    is_starred = lambda x: x >= args.starred
    blast_results['starred'] = blast_results['pident'].apply(is_starred)

    # merge blast results with seq_info - do this early so that
    # refseqs not represented in the blast results are discarded in
    # the merge.
    blast_results = pd.merge(
        # blast_results[['qseqid', 'sseqid', 'pident', 'starred', 'coverage']],
        blast_results[['qseqid', 'sseqid', 'pident', 'starred']],
        seq_info[['tax_id', 'accession']],
        how='left', left_on='sseqid', sort=False, right_index=True)

    # merge with taxonomy to associate original tax_ids with tax_ids
    # at the specified rank.
    blast_results = pd.merge(blast_results,
                             taxonomy[[args.rank]],
                             how = 'left',
                             left_on = 'tax_id',
                             right_index = True)

    # merge with taxonomy again to get tax_names
    blast_results = pd.merge(blast_results,
                             taxonomy[['tax_name']],
                             how = 'left',
                             left_on = args.rank,
                             right_index = True)

    # output to details.csv.bz2
    if args.out_details:
        if args.details_full:
            with args.out_details as out_details:
                blast_results.to_csv(out_details,
                                     header=True,
                                     index=False)
        else:
            print 'details summary not implemented yet'
            return

    ## concludes our blast results details, now for some grouping and summarizing

    # Create a Series containing each of sequence names (ie, qseqids) with
    # values indicating no blast hits (0), < min_identity (1), or >=
    # max_idenity (2).
    no_hit_idx = blast_results.sseqid.isnull()

    not_assigned = blast_results[no_hit_idx]
    not_assigned['assignment'] = '[no blast_result]'

    tax_dict = {i:t.to_dict() for i,t in taxonomy.fillna('').iterrows()}

    def assign(df):
        ids_stars = df.groupby(by=['tax_id', 'starred']).groups.keys()
        condensed = sequtils.condense_assignment(ids_stars, tax_dict)
        assignment = sequtils.compound_assignment(condensed, tax_dict)
        df['assignment'] = assignment
        return df

    assigned = blast_results[~ no_hit_idx].groupby(['qseqid']).apply(assign)

    assigned = pd.concat([assigned, not_assigned])

    if args.weights:
        weights = read_csv(args.weights, header=None, names=['qseqid', 'weight'])
        assigned = pd.merge(assigned, weights, how='left', on='qseqid')
    else:
        assigned['weight'] = 1

    # group by assignment
    counts = assigned.groupby(['assignment']).weight.sum()
    counts.name = 'reads'

    output = pd.DataFrame(counts)
    output['freq'] = output / output.sum()
    output['clusters'] = assigned.groupby(['assignment']).size()

    output.sort('reads', ascending=False, inplace=True)

#    output.to_csv(args.outfile)

    # next steps...
    # - create names for assignments
    # - calculate and apply copy number correction factor
    # - add assignment_ids
    # - assemble details


