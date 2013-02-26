"""Run-length encode a fasta file"""

import logging
import csv
import sys

from itertools import imap, chain

from bioy_pkg.sequtils import homoencode, to_ascii, fastalite
from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('infiles',
            nargs = '*',
            default = [sys.stdin],
            type = Opener(),
            help = 'Input fasta file')
    parser.add_argument('-o','--outfile',
            type = Opener('w'),
            default = sys.stdout,
            help = 'Name of output file; default: %(default)s')
    parser.add_argument('-r','--rlefile',
            type = Opener('w'),
            help = """Name of output file for run length encoding; default is to
                      append .csv.bz2 to --outfile basename.""")

def action(args):
    counts = []
    for seq in chain.from_iterable(imap(fastalite, args.infiles)):
        seqstr, count = homoencode(seq.seq)
        assert len(seqstr) == len(count)
        args.outfile.write('>{}\n{}\n'.format(seq.description, seqstr))
        counts.append([seq.description, to_ascii(count)])

    if args.rlefile:
        writer = csv.writer(args.rlefile)
        writer.writerow(['name', 'rle'])
        writer.writerows(counts)

