"""
Classify sequences by grouping blast output by matching taxonomic names

Optional grouping by specimen and query sequences
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
        help='tax table of taxids and species names')

    parser.add_argument('--all-one-group',
                        dest = 'all_one_group',
            action = 'store_true', default=False,
            help = """If --map is not provided, the default behavior is to treat
                    all reads as one group; use this option to treat
                    each read as a separate group [%(default)s]""")
    parser.add_argument(
        '--starred', default = 100.0, metavar='PERCENT', type = float,
        help = """Names of organisms for which at least one reference
        sequence has pairwise identity with a query sequence of at
        least PERCENT will be marked with an asterisk [%(default)s]""")
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
    parser.add_argument('--rank',
            metavar='RANK',
            help = 'Rank at which to classify. Default: "%(default)s"',
            default = 'species')
    ### csv.Sniffer.has_header is *not* reliable enough
    parser.add_argument('--has-header', action = 'store_true',
            help = 'specify this if blast data has a header')


class Assignment(object):
    def __init__(self, colnames):
        self.colnames = colnames
        self.assignments = {}

    def assign(self, df):
        assignment = frozenset(df.groupby(self.colnames).groups.keys())
        h = hash(assignment)
        self.assignments[h] = assignment
        return h


def read_csv(filename, compression=None, **kwargs):
    """read a csv file using pandas.read_csv with compression defined by
    the file suffix unless provided.

    """

    suffixes = {'.bz2': 'bz2', '.gz': 'gzip'}
    kwargs['compression'] = compression or suffixes.get(path.splitext(filename)[-1])

    return pd.read_csv(filename, **kwargs)


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
    blast_results = read_csv(args.blast_file, header=None)
    fieldnames = None if args.has_header else sequtils.BLAST_HEADER
    blast_results.rename(columns=dict(zip(blast_results.columns, fieldnames)),
                         inplace=True)

    # get rows corresponding to sequences with no blast hits
    no_hit_rows = blast_results['sseqid'].isnull()
    no_hits = blast_results[no_hit_rows]

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

    assignments = grouped_by_query.apply(assigner.assign)

    print assignments

    # for k, v in assigner.assignments.items():
    #     print sorted(v)

    # next steps...
    # - concatenate with no_hits (corresponding to an assignment of 0? nan?)
    # - create a DataFrame indexed by qseqid with assignments and additional column 'weights'
    # - group by assignment and get sum of weights
    # - create names for assignments
    # - calculate and apply copy number correction factor
    # - order assignments by mass desc and add assignment_id
    # - assemble details


