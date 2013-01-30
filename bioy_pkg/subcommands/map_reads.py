"""
Map fasta reads to csv sample file
"""

import sys
import logging
import re
import csv

from bioy_pkg.utils import Opener, opener
from bioy_pkg.sequtils import fastalite

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('infiles',
            nargs = '+',
            help = 'Input fasta files')
    parser.add_argument('-b', '--barcode',
            default = '/(.*)\.',
            help = 'Regex of file name to barcode label. Default: %(default)s')
    parser.add_argument('-t', '--type',
            choices = ['ion', '454'],
            default = 'ion',
            help = 'Read type. Default: %(default)s')
    parser.add_argument('-i', '--id',
            default = '1',
            help = 'Plate number')
    parser.add_argument('-o','--out',
            default = sys.stdout,
            type = Opener('w'),
            help = 'csv(.bz2) output file')

def action(args):
    pzr = 'i' if args.type == 'ion' else 'c'
    pzr += args.id
    writer = csv.writer(args.out)
    for infile in args.infiles:
        sample = '{}{}'.format(
                pzr, ''.join(re.search(args.barcode, infile).groups()))
        for fasta in fastalite(opener(infile)):
            writer.writerow([fasta.id, sample])

