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
                        default=sys.stdin,
                        nargs='?',
                        type=Opener(),
                        help='input csv file')
    parser.add_argument('-o', '--outfile',
                        help='HDF5 file [use basename of input file by default]')
    parser.add_argument('-d', '--outdir',
                        help="""Optional output directory. Ignored
                        if -o/--outfile is specified.""")
    parser.add_argument('--fieldnames',
                        help='comma-delimited list of field names.')
    parser.add_argument('-H', '--no-header', action='store_true', default=False,
                        help="""indicate that the input file has no
                        header row. Uses value of --fieldnames if provided.""")
    parser.add_argument('-k', '--key', default='data',
                        help="""A label identifing this table in the
                        data store "[%(default)s]" """)
    parser.add_argument('-c', '--no-compress', action='store_false', default=True,
                        dest='compress', help="""Don't compress data store.""")


def action(args):

    try:
        import tables
    except ImportError, e:
        sys.exit('This subcommand requires PyTables and dependencies '
                 '(see README and requirements.txt)')

    if args.outfile:
        outfile = args.outfile
    else:
        indir, infile = path.split(args.infile.name)
        outdir = args.outdir or indir
        outfile = path.join(outdir, path.splitext(infile)[0]) + '.hdf5'

    tab = pd.read_csv(args.infile,
                      header=None if args.no_header else 0)

    fieldnames = args.fieldnames.split(',') if args.fieldnames else None
    if fieldnames:
        tab.rename(columns=dict(zip(tab.columns, fieldnames)), inplace=True)

    # tab[[2, 6]] = tab[[2, 6]].astype(np.float16)
    # tab[[3, 4, 5]] = tab[[3, 4, 5]].astype(np.int16)
    # tab[[0, 1]] = tab[[0, 1]].astype(pd.Categorical)

    if args.verbosity > 1:
        print tab.head()
        print tab.dtypes

    hdfargs = {}
    if args.compress:
        hdfargs.update(dict(complevel=9, complib='bzip2'))

    tab.to_hdf(outfile, args.key, **hdfargs)
