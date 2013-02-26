"""Run-length encode a fasta file"""

import logging
import sys

from itertools import imap, chain
from csv import DictWriter

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
            type = lambda f: DictWriter(Opener('w')(f), fieldnames=['name', 'rle']),
            help = """Name of output file for run length encoding; default is to
                      append .csv.bz2 to --outfile basename.""")

def action(args):
    for seq in chain.from_iterable(imap(fastalite, args.infiles)):
        seqstr, count = homoencode(seq.seq)
        assert len(seqstr) == len(count)
        args.outfile.write('>{}\n{}\n'.format(seq.description, seqstr))
        if args.rlefile:
            args.rlefile.writerow(dict({'name':seq.description, 'rle':to_ascii(count)}))

