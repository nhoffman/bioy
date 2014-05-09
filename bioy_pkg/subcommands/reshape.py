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
