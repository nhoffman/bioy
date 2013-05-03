"""
Calculate all Smith-Waterman pairwise distances among sequences.
"""

import logging
import sys
import pprint
import csv

from itertools import islice, chain, groupby, imap
from operator import itemgetter

from bioy_pkg.sequtils import fastalite, all_pairwise
from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('fasta',
                        default=sys.stdin,
                        type=Opener('r'),
                        help='sequences in fasta format')
    parser.add_argument('-o', '--out',
                        default=sys.stdout,
                        type=Opener('w'),
                        help='output file (default stdout)')
    parser.add_argument('-d', '--distance', action='store_true',
                        default=False, help = 'Calculate distance rather than identity.')

def action(args):

    seqs = fastalite(args.fasta)
    pairs = all_pairwise(seqs)

    writer = csv.writer(args.out)
    writer.writerow(['query', 'target', 'identity'])

    if args.distance:
        pairs = ((q, t, 1 - i) for q, t, i in pairs)

    writer.writerows(pairs)
