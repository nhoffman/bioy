"""
Turn a fasta file into a csv
"""

import csv
import logging
import operator
import sys

from bioy_pkg.sequtils import fastalite
from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('fasta',
            nargs = '?',
            default = sys.stdin,
            metavar = 'FILE',
            type = Opener(),
            help = 'A fasta file')
    parser.add_argument('--get',
            action = 'append',
            help = 'columname[:newname]')
    parser.add_argument('--out',
            metavar = 'FILE',
            type = Opener('w'),
            default = sys.stdout,
            help = 'output csv file columns: [id,description,seq]')

def action(args):
    fieldnames = args.get or ['id','description','seq']
    # make into [[columname,newname] ...]
    fieldnames = [f.split(':') for f in fieldnames]
    fieldnames = [f * (2 if len(f) == 1 else 1) for f in fieldnames]
    out = csv.DictWriter(args.out,
                         fieldnames = map(operator.itemgetter(1), fieldnames),
                         extrasaction = 'ignore')
    out.writeheader()
    for f in fastalite(args.fasta):
        f = f._asdict()
        out.writerow({v:f[k] for k,v in fieldnames})

