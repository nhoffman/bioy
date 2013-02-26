"""Run-length decode a fasta file"""

import logging
import sys

from bioy_pkg.sequtils import homodecode, from_ascii, fastalite
from bioy_pkg.utils import Opener, Csv2Dict

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('infile',
            type = Opener(),
            help = 'Input fasta file')
    parser.add_argument('rlefile',
            type = Csv2Dict(),
            help = 'csv file (may be bzip encoded) containing columns "name","rle"')
    parser.add_argument('-o','--outfile',
            type = Opener('w'),
            default = sys.stdout,
            help = 'Name of output file')

def action(args):
    for s in fastalite(args.infile):
        rle = from_ascii(args.rlefile[s.description])
        assert len(s.seq) == len(rle)
        args.outfile.write('>{}\n{}\n'.format(s.description, homodecode(s.seq, rle)))
