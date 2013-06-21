"""reverse complement rle and non-rle sequences"""

import logging
import sys

from csv import DictWriter

from bioy_pkg.sequtils import fastalite
from bioy_pkg.utils import Opener, Csv2Dict

log = logging.getLogger(__name__)

rle_fieldnames = ['name', 'rle']

def build_parser(parser):
    parser.add_argument('infile',
            type = Opener(),
            help = 'Input fasta file')
    parser.add_argument('rlefile',
            nargs = '?',
            type = Csv2Dict(value = 'rle',
                            fieldnames = rle_fieldnames),
            help = 'csv file (may be bzip encoded) containing columns "name","rle"')
    parser.add_argument('-O', '--out-rle',
            type = lambda f: DictWriter(Opener('w')(f), fieldnames = rle_fieldnames),
            help = 'reversed rlefile')
    parser.add_argument('-o','--out-fasta',
            type = Opener('w'),
            default = sys.stdout,
            help = 'Name of output file')

def action(args):
    seqs = fastalite(args.infile)

    for s in seqs:
        seq = ''.join(reversed(s.seq))
        args.out_fasta.write('>{}\n{}\n'.format(s.description, seq))

    if args.rlefile and args.out_rle:
        for n,r in args.rlefile.items():
            rle = ''.join(reversed(r))
            args.out_rle.writerow(dict(name = n, rle = rle))
