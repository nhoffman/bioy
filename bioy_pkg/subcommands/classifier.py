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

import pandas as pd

from bioy_pkg import sequtils
from bioy_pkg.utils import Opener, opener, Csv2Dict, groupbyl

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('blast_file',
            help = 'CSV tabular blast file of query and subject hits.')
    parser.add_argument('-s', '--seq-info',
            required = True,
            metavar = 'CSV',
            help = 'File mapping reference seq name to tax_id')
    parser.add_argument('-t', '--taxonomy',
            required = True,
            metavar = 'CSV',
            help = 'tax table of taxids and species names')

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
    parser.add_argument('-w', '--weights',
            metavar = 'CSV',
            type = Opener(),
            help = 'columns: name, weight')
    ### csv.Sniffer.has_header is *not* reliable enough
    parser.add_argument('--has-header', action = 'store_true',
            help = 'specify this if blast data has a header')


def action(args):
    # format blast data and add additional available information
    fieldnames = None if args.has_header else sequtils.BLAST_HEADER

    seq_info = pd.read_csv(args.seq_info,
                           dtype={'seqname': str,
                                  'tax_id': str,
                                  'accession': str,
                                  'description': str,
                                  'length': int,
                                  'ambig_count': int,
                                  'is_type': str,
                                  'rdp_lineage': str})
    seq_info.rename(columns=dict(seqname='sseqid'), inplace=True)

    taxonomy = pd.read_csv(args.taxonomy, dtype=str)

    compression = dict(bz2='bz2', gz='gzip').get(path.splitext(args.blast_file)[-1])
    blast_results = pd.read_csv(args.blast_file, header=None, compression='bz2')
    blast_results.rename(columns=dict(zip(blast_results.columns, fieldnames)),
                         inplace=True)

    # merge blast results with seq_info - do this first so that
    # refseqs not represented in the blast results are discarded.
    hits = pd.merge(blast_results[['qseqid', 'sseqid', 'pident', 'coverage']],
                    seq_info[['sseqid', 'tax_id', 'accession']],
                    how='left', on='sseqid', sort=False)

    # merge with taxonomy to map original tax_ids with tax_ids at the
    # specified rank.
    hits = pd.merge(hits, taxonomy[['tax_id', args.rank]], how='left',
                    on='tax_id', sort=False)

    # merge with taxonomy again to get tax_names
    hits = pd.merge(hits, taxonomy[['tax_id', 'tax_name']], how='left',
                    left_on=args.rank, right_on='tax_id')
    del hits['tax_id_y']

    groups = hits.groupby(['qseqid'])

    for qseqid, group in groups:
        print group.tax_name.unique()

    # what are nan results?


