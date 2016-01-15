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
Turn a csv file into a fasta file specifying two columns
"""

import logging
import sys

from csv import DictReader

from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('csv',
            nargs = '?',
            default = sys.stdin,
            metavar = 'FILE',
            type = Opener(),
            help = 'A csv file with at least one column')
    parser.add_argument('--columns',
            default = '1,2',
            help = 'Comma-delimited list of column names or numbers')
    parser.add_argument('--out',
            metavar = 'FILE',
            type = Opener('w'),
            default = sys.stdout,
            help = 'output fasta file')

def action(args):
    cols = [c.strip() for c in args.columns.split(',')]

    length = len(cols)

    reader = DictReader(args.csv)

    if all(c.isdigit() for c in cols):
        name = reader.fieldnames[int(cols[0 % length])]
        seq = reader.fieldnames[int(cols[1 % length])]
    else:
        name = cols[0]
        seq = cols[1]

    for row in reader:
        args.out.write('>{}\n{}\n'.format(row[name], row[seq]))

