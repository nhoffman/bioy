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
convert a tsv file to a csv with an optional split/add columns feature
"""

import logging
import numpy
import pandas
import sys

from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('csv',
            default = sys.stdin,
            nargs = '?',
            type = Opener(),
            help = 'input tsv file')
    parser.add_argument('-o', '--out',
            type = Opener('w'),
            default = '/dev/stdout',
            help = 'csv file')
    parser.add_argument('--value',
            help = 'value to pivot on for each row and column')
    parser.add_argument('--rows',
            help = 'comma delimited list of values')
    parser.add_argument('--cols',
            help = 'comma delimited list of values')
    parser.add_argument('--fill-value',
            help = 'if no pivot value')

def action(args):
    df = pandas.read_csv(args.csv)
    df = pandas.pivot_table(df, values=args.value,
                                rows=args.rows.split(','),
                                cols=args.cols.split(','),
                                fill_value=args.fill_value,
                                aggfunc=numpy.sum)
    df.to_csv(args.out)
