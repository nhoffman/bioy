"""
Create seq info file from fasta file
"""

import logging
import sys

from csv import DictWriter

from bioy_pkg.sequtils import fastalite, count_ambiguous
from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('fasta',
            type = Opener(),
            nargs = '?',
            default = sys.stdin,
            help = 'input file containing raw reads')
    parser.add_argument('--fieldnames',
            type = lambda f: [n.strip() for n in f.split(',')],
            default = ['seqname', 'description', 'ambig_count', 'length'],
            help = 'seq info comma seperated field names.  Default = %(default)s')
    parser.add_argument('-o', '--out',
            type = Opener('w'),
            default = sys.stdout,
            help = 'fasta output file')

def action(args):
    fasta = fastalite(args.fasta)

    info = DictWriter(args.out, fieldnames = args.fieldnames, extrasaction = 'ignore')

    info.writeheader()

    for f in fasta:
        info.writerow({
            'seqname':f.id,
            'description':f.description,
            'ambig_count':count_ambiguous(f.seq),
            'length':len(f.seq)
            })
