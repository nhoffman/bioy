# This file is part of Bioy
#
#    Bioy is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Bioy is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Bioy.  If not, see <http://www.gnu.org/licenses/>.

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

