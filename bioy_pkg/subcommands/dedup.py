"""
Deduplicate sequences by coalescing identical substrings

Use seqmagick convert --deduplicate sequences for fast deduplication
of identical sequences.
"""

import logging
import sys
from itertools import groupby

from bioy_pkg.deduplicate import dedup
from bioy_pkg.sequtils import fastalite
from bioy_pkg.utils import Csv2Dict, opener, Opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('sequences',
            type = lambda f: list(fastalite(opener(f))),
            help = 'input fasta file')
    parser.add_argument('-i', '--seq-info',
            type = Csv2Dict('seqname'),
            help = """csv file containing column "seqname" plus
                      another column for grouping sequences prior to deduplication""")
    parser.add_argument('--primary-group',
            help='string specifying column in seq_info to use for grouping [default %(default)s]',
            default = 'species')
    parser.add_argument('--secondary-group',
            help = """string specifying column in seq_info to use for grouping
                      if primary_group is undefined for a given row [default %(default)s]""",
            default = 'tax_id')
    parser.add_argument('-o','--out',
            type = Opener('w'),
            default = sys.stdout,
            help = 'Output fasta file.')

def action(args):
    def groupfun(seq):
        lineage = args.seq_info[seq.description]
        if args.secondary_group:
            return lineage[args.primary_group] or lineage[args.secondary_group]
        else:
            return lineage[args.primary_group]

    if args.seq_info:
        for grp, seqs in groupby(sorted(args.sequences, key = groupfun), groupfun):
            seqs = list(seqs)
            deduped = dedup([s.seq for s in seqs])
            for i in deduped.keys():
                args.out.write('>{}\n{}\n'.format(seqs[i].description, seqs[i].seq))
