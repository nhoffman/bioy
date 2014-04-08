"""
Calculate all Smith-Waterman pairwise distances among sequences.
"""

import logging
import sys
import csv

from itertools import groupby
from numpy import median
from operator import itemgetter

from bioy_pkg.sequtils import fastalite, all_pairwise
from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('fasta', default=sys.stdin, type=Opener('r'),
                        help='sequences in fasta format')
    parser.add_argument('-o', '--out', default=sys.stdout, type=Opener('w'),
                        help='output file (default stdout)')
    parser.add_argument('--median-out', type = Opener('w'),
                        help = """median score of pairwise alignments""")
    parser.add_argument('-d', '--distance', action='store_true', default=False,
                        help = 'Calculate distance rather than identity.')

def action(args):

    seqs = fastalite(args.fasta)
    pairs = list(all_pairwise(seqs))

    if args.distance:
        pairs = [(q, t, 1 - i) for q, t, i in pairs]

    writer = csv.writer(args.out)
    writer.writerow(['query', 'target', 'identity'])
    writer.writerows(pairs)

    if args.median_out:
        pairs += map(itemgetter(1,0,2), pairs)
        pairs = sorted(pairs, key = itemgetter(0))
        pairs = groupby(pairs, key = itemgetter(0))

        writer = csv.writer(args.median_out)
        writer.writerow(['query', 'median'])

        for seqname, group in pairs:
            scores = map(itemgetter(2), group)
            med = median(scores)
            log.info('{}: median({}) = {}'.format(seqname, scores, med))
            writer.writerow([seqname, med])

