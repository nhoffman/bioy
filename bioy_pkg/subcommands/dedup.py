"""
Deduplicate sequences by coalescing identical substrings

Use seqmagick convert --deduplicate sequences for fast deduplication
of identical sequences.
"""

import logging
import sys

from csv import DictWriter
from itertools import groupby, chain

from bioy_pkg.deduplicate import dedup
from bioy_pkg.sequtils import fastalite
from bioy_pkg.utils import Csv2Dict, opener, Opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('sequences',
            type = lambda f: fastalite(opener(f), readfile = False),
            help = 'input fasta file')
    parser.add_argument('-i', '--seq-info',
            type = Csv2Dict('seqname'),
            help = """csv file containing column "seqname" plus
                      another column for grouping sequences prior to deduplication""")
    parser.add_argument('--primary-group',
            help = 'string specifying column in seq_info to use for grouping [default %(default)s]',
            default = 'species')
    parser.add_argument('--secondary-group',
            help = """string specifying column in seq_info to use for grouping
                      if primary_group is undefined for a given row [default %(default)s]""",
            default = 'tax_id')
    parser.add_argument('-O', '--out-info',
            type = Opener('w'),
            help = 'deduplicate seq info file')
    parser.add_argument('-o','--out',
            type = Opener('w'),
            default = sys.stdout,
            help = 'Output fasta file')

def dedups(seqs):
    d = dedup([s.seq for s in seqs]).keys()
    d = map(lambda i: seqs[i], d)
    return d

def action(args):
    def seq_group(seq):
        if args.seq_info:
            lineage = args.seq_info[seq.id]
            if args.secondary_group:
                return lineage[args.primary_group] or lineage[args.secondary_group]
            else:
                return lineage[args.primary_group]
        else:
            return 'all'

    seqs = sorted(args.sequences, key = seq_group)
    seqs = groupby(seqs, seq_group)
    seqs = (dedups(list(s)) for _,s in seqs)
    seqs = chain(*seqs)

    if args.out_info and args.seq_info:
        fieldnames = next(iter(args.seq_info.values())).keys()
        info_out = DictWriter(args.out_info, fieldnames = fieldnames)
        info_out.writeheader()

    for s in seqs:
        args.out.write('>{}\n{}\n'.format(s.description, s.seq))

        if args.out_info and args.seq_info:
            info_out.writerow([args.seq_info[s.id]])

