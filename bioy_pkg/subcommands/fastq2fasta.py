"""Convert fastq file to fasta plus quality scores (Not completed yet)"""

import logging
import sys

from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('infile',
            nargs = '?',
            default = sys.stdin,
            type = Opener(),
            help = 'Input fastq file')
    parser.add_argument('-o','--outfile',
            type = Opener('w'),
            help = 'Name of output file; default is to append .fasta to basename.')
    parser.add_argument('-q','--quality-scores',
            type = Opener('w'),
            help = 'Name of output file for quality scores; default is to append .qual to basename.')

def action(args):
    pass
