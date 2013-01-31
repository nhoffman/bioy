"""Run-length encode a fasta file"""

import logging
import csv
import sys

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
            help = 'Name of output file; default is to append _rle.fasta to basename.')
    parser.add_argument('-r','--rlefile',
            type = Opener('w'),
            help = 'Name of output file for run length encoding; default is to append _rle.csv.bz2 to basename.')

def action(args):
    seqs = fastalite(args.infile)

    outfile = args.outfile or opener(args.infile.replace('.fasta', '_rle.fasta', 'w'))
    rlefile = args.rlefile or opener(args.infile.replace('.fasta', '_rle.csv.bz2'), 'w')

    writer = csv.writer(rlefile)
    writer.writerow(['name','rle'])
    for seq in seqs:
        seqstr, counts = homoencode(seq.seq)
        assert len(seqstr) == len(counts)
        outfile.write('>{}\n{}\n'.format(seq.description, seqstr))
        writer.writerow([seq.description, to_ascii(counts)])
