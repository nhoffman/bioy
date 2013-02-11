"""Run-length encode a fasta file"""

import logging
import csv
import sys
from concurrent.futures import ProcessPoolExecutor

from os.path import splitext

from bioy_pkg.sequtils import homoencode, to_ascii, fastalite
from bioy_pkg.utils import Opener, opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('infile',
            nargs = '?',
            default = sys.stdin,
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
    if args.rlefile:
        rlefile = args.rlefile
    elif args.outfile and args.outfile is not sys.stdout:
        rlefile = opener('{}.{}'.format(splitext(args.outfile.name)[0], 'csv.bz2'), 'w')
    else:
        rlefile = None

    counts = []
    for seq in fastalite(args.infile):
        seqstr, count = homoencode(seq.seq)
        assert len(seqstr) == len(count)
        args.outfile.write('>{}\n{}\n'.format(seq.description, seqstr))
        counts.append([seq.description, to_ascii(count)])

    if rlefile:
        writer = csv.writer(rlefile)
        writer.writerow(['name', 'rle'])
        writer.writerows(counts)

