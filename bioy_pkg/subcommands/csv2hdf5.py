"""
Convert a csv file to HDF5
"""

from os import path
import csv
import logging
import re
import sys

import pandas as pd

from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('infile',
            default = sys.stdin,
            nargs = '?',
            type = Opener(),
            help = 'input csv file')
    parser.add_argument('-o', '--outfile',
                        help = 'HDF5 file [use basename of input file by default]')
    parser.add_argument('-d', '--outdir',
                        help = 'Optional output directory. Ignored if -o/--outfile is specified.')
    parser.add_argument('--fieldnames', help='comma-delimited list of field names.')
    parser.add_argument('--headless', action='store_true', default=False,
                        help="""indicate that the input file has no
                        header row. Uses value of --fieldnames if provided.""")
    parser.add_argument('-k', '--key', default='data',
                        help='A label identifing this table in the data store "[%(default)s]"')
    parser.add_argument('-c', '--compress', action='store_true', default=False,
                        help="""Compress data store.""")
    parser.add_argument('-n', '--dry-run', action='store_true', default=False,
                        help="""Print a representation of the input data and exit""")


def action(args):


    try:
        import tables
    except ImportError, e:
        sys.exit('This subcommand requires PyTables and dependencies (see README and requirements.txt)')

    if args.outfile:
        outfile = args.outfile
    else:
        indir, infile = path.split(args.infile.name)
        outdir = args.outdir or indir
        outfile = path.join(outdir, path.splitext(infile)[0]) + '.hdf5'

    tab = pd.read_csv(args.infile,
                      header=None if args.headless else 0)

    fieldnames = args.fieldnames.split(',') if args.fieldnames else None
    if fieldnames:
        tab.rename(columns=dict(zip(tab.columns, fieldnames)), inplace=True)

    if args.dry_run:
        print tab.head()
        print tab.dtypes
        sys.exit(0)

    tab.to_hdf(outfile, args.key)
