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
from csv import DictReader, DictWriter
from collections import defaultdict
from itertools import groupby
from math import ceil
from operator import itemgetter

import numpy as np
import pandas as pd

from bioy_pkg import sequtils
from bioy_pkg.utils import Opener, opener, Csv2Dict, groupbyl

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
        '-o', '--out', default=sys.stdout, type=Opener('w'), metavar='FILE',
        help="Classification results.")
    parser.add_argument(
        '-d', '--details', type=Opener('w'), metavar='FILE',
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
    parser.add_argument(
        '--group-label', metavar = 'LABEL', default = 'all',
        help = """Single group label for reads""")

    # options for details file
    parser.add_argument(
        '--details-identity', metavar='PERCENT', type=float, default=90.0,
        help='Minimum identity to include blast hits in details file')
    parser.add_argument(
        '--details-full', action='store_true', default='False',
        help='do not limit out_details to only larget cluster per assignment')


class Assignment(object):
    """An instance of the Assignment class is used to create a dictionary
    of group keys. Each time the assign() method is called, the
    provided DataFrame 'df' is grouped according to self.colnames and
    the group keys are added to self.assignment keyed by a hash.

    """

    def __init__(self, colnames):
        """'colnames' is a list of column names on which to group the input of
        the 'assign' method.

        """
        self.colnames = colnames
        self.assignments = {}

    def assign(self, df):
        """Store the group keys resulting from grouping 'df' by self.colnames
        in self.assignments and return the hash.

        """
        assignment = frozenset(df.groupby(self.colnames).groups.keys())
        assignment_key = hash(assignment)
        self.assignments[assignment_key] = assignment
        return assignment_key


def assign_names(assignment, no_result='[no blast result]'):
    """Mockup of a function to convert a frozenset containing (tax_id,
    is_starred) tuples into a taxonomic name.

    """

    if assignment:
        try:
            return '-'.join(sorted({str(tax_id) for tax_id, _ in assignment}))
        except TypeError:
            print assignment
            sys.exit()
    else:
        return no_result


def read_csv(filename, compression=None, **kwargs):
    """read a csv file using pandas.read_csv with compression defined by
    the file suffix unless provided.

    """

    suffixes = {'.bz2': 'bz2', '.gz': 'gzip'}
    kwargs['compression'] = compression or suffixes.get(path.splitext(filename)[-1])

    return pd.read_csv(filename, **kwargs)


# details columns (maybe)
# ['specimen', 'assignment', 'assignment_id', 'qseqid', 'sseqid',
# 'pident', 'coverage', 'ambig_count', 'accession', 'tax_id',
# 'tax_name', 'target_rank', 'rank', 'hi', 'low']

def action(args):

    seq_info = read_csv(args.seq_info,
                        dtype={'seqname': str,
                               'tax_id': str,
                               'accession': str,
                               'description': str,
                               'length': int,
                               'ambig_count': int,
                               'is_type': str,
                               'rdp_lineage': str})

    seq_info.rename(columns=dict(seqname='sseqid'), inplace=True)

    taxonomy = read_csv(args.taxonomy, dtype=str)

    # format blast data and add additional available information
    blast_results = read_csv(
        args.blast_file, header=0 if args.has_header else None)

    if not args.has_header:
        # assert len(blast_results.columns) == len(sequtils.BLAST_HEADER)
        blast_results.rename(
            columns=dict(zip(blast_results.columns, sequtils.BLAST_HEADER)),
            inplace=True)

    # debug: keep track of original qseqids here
    qseqids_in = {x for x in blast_results.qseqid}

    # Create a Series containing sequence names (ie, qseqids) with no
    # blast hits; represent this state with a value of 0.
    no_hit_rows = blast_results['sseqid'].isnull()
    no_hits = pd.Series(0, index=blast_results.qseqid[no_hit_rows])

    # ... and remove them from blast_results
    blast_results = blast_results[no_hit_rows.apply(lambda x: not x)]

    # add a column indicating if a hit meets the threshold for
    # starring.
    blast_results['starred'] = blast_results['pident'].apply(
        lambda x: x >= args.starred)

    # merge blast results with seq_info - do this first so that
    # refseqs not represented in the blast results are discarded.
    hits = pd.merge(
        # blast_results[['qseqid', 'sseqid', 'pident', 'starred', 'coverage']],
        blast_results[['qseqid', 'sseqid', 'pident', 'starred']],
        seq_info[['sseqid', 'tax_id', 'accession']],
        how='left', on='sseqid', sort=False)

    # merge with taxonomy to associate original tax_ids with tax_ids
    # at the specified rank.
    hits = pd.merge(hits, taxonomy[['tax_id', args.rank]], how='left',
                    on='tax_id', sort=False)

    # merge with taxonomy again to get tax_names
    hits = pd.merge(hits, taxonomy[['tax_id', 'tax_name']], how='left',
                    left_on=args.rank, right_on='tax_id')

    # rename the column containing tax_id's at the specified rank to
    # 'tax_id' and delete the duplicate column 'tax_id_y'; tax_id_x
    # now refers to the original tax_ids for the reference sequences.
    hits.rename(columns={args.rank: 'tax_id'}, inplace=True)
    del hits['tax_id_y']

    # stores a hash of each frozenset containing {(tax_name, starred), ...}
    assigner = Assignment(['tax_id', 'starred'])

    grouped_by_query = hits.groupby(['qseqid'])

    # Determine assignment names and concatenate with sequences having
    # no blast hit.
    assigned = pd.concat([no_hits, grouped_by_query.apply(assigner.assign)])
    assigned.name = 'assignment_key'

    assert set(assigned.index) == qseqids_in

    assignments = pd.DataFrame(assigned)
    assignments['qseqid'] = assignments.index

    if args.weights:
        weights = read_csv(args.weights, header=None, names=['qseqid', 'weight'])
        assignments = pd.merge(assignments, weights, how='left', on='qseqid')
    else:
        assignments['weight'] = 1

    # add assignments
    assignments['assignment'] = assignments.assignment_key.map(
        lambda x: assign_names(assigner.assignments.get(x)))

    # group by assignment
    counts = assignments.groupby(['assignment']).weight.sum()
    counts.name = 'reads'

    output = pd.DataFrame(counts)
    output['freq'] = output/output.sum()
    output['clusters'] = assignments.groupby(['assignment']).size()

    output.sort('reads', ascending=False, inplace=True)

    output.to_csv(args.out)

    # next steps...
    # - create names for assignments
    # - calculate and apply copy number correction factor
    # - add assignment_ids
    # - assemble details


