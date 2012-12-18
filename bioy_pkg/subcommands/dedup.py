"""
Deduplicate sequences by coalescing identical substrings

Use seqmagick convert --deduplicate sequences for fast deduplication
of identical sequences.
"""

import logging
import sys
import argparse
import csv
from itertools import groupby, ifilter, islice
from random import shuffle
from collections import defaultdict

from bioy_pkg.deduplicate import dedup
from bioy_pkg.sequtils import fastalite
from bioy_pkg.utils import csv2dict

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('infile', type = argparse.FileType('r'), help = 'input fasta file')
    parser.add_argument('-i','--seq-info', type = csv2dict('seqname'),
                        help='csv file containing column "seqname" plus another column for grouping sequences prior to deduplication')
    parser.add_argument('--primary-group', help='string specifying column in seq_info to use for grouping [default %(default)s]',
                        default = 'species')
    parser.add_argument('--secondary-group', help='string specifying column in seq_info to use for grouping if primary_group is undefined for a given row [default %(default)s]',
                        default = 'tax_id')
    parser.add_argument('-o','--outfile', type=argparse.FileType('w'), default = sys.stdout,
                        help='Output fasta file.')

def action(args):
    seqs = list(fastalite(args.infile))

    def groupfun(seq):
        lineage = args.seq_info[seq.id]
        if args.secondary_group:
            return lineage[args.primary_group] or lineage[args.secondary_group]
        else:
            return lineage[args.primary_group]

    if args.seq_info:
        for grp, seqs in groupby(sorted(seqs, key = groupfun), groupfun):
            seqs = list(seqs)
            deduped = dedup([s.seq for s in seqs])
            for i in deduped.keys():
                args.outfile.write('>{}\n{}\n'.format(seqs[i].id, seqs[i].seq))
