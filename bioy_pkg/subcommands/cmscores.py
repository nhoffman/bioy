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
Convert raw cmalign alignment scores to csv format.

Discard the first column ("seq idx") as well as the "elapsed" column
if present.
"""

import logging
import sys
import csv

from itertools import imap, groupby, ifilter

from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)

CMALIGN_SCORE_FIELDS = ['name','len','total','struct','prob']

def build_parser(parser):
    parser.add_argument('infile',
            type = Opener('r'),
            default = sys.stdin,
            nargs = '?',
            help = 'File containing stdout of call to cmalign')
    parser.add_argument('-o','--outfile',
            type = Opener('w'),
            default = sys.stdout,
            help = 'Output file in csv format (default is stdout)')

def read_scores(infile):

    """
    Return an iterator of tuples in open file object ``infile``. Example input file:

    # seq idx  seq name          len     total    struct  avg prob      elapsed
    # -------  --------------  -----  --------  --------  --------  -----------
            1  FUM0LCO01B97N3    248    385.89     43.47     1.000  00:00:00.07
            2  FUM0LCO01EZTWO    248    385.89     43.47     1.000  00:00:00.07
            3  FUM0LCO01ALN4O    241    374.25     44.48     0.990  00:00:00.07
            4  FUM0LCO01D7IKL    242    385.80     45.69     1.000  00:00:00.08

    Note that elapsed seems to be omitted when cmalign is run with mpi
    """
    groups = groupby(infile, lambda line: line.startswith('#'))

    is_comment, data = next(groups)

    assert is_comment is True

    # parse header from comments
    header, = [line.strip() for line in data if line.startswith('# seq idx')]

    # Throw out the first column; also throw out the last column if
    # it's called "elapsed" (only present when cmalign is run without
    # mpi).
    if header.endswith('elapsed'):
        rowmap = lambda row: row.split()[1:-1]
    else:
        rowmap = lambda row: row.split()[1:]

    # parse data
    is_comment, data = next(groups)

    assert is_comment is False

    data = ifilter(lambda row: bool(row.strip()), data)

    data = imap(rowmap, data)

    return data

def action(args):
    writer = csv.writer(args.outfile)
    writer.writerow(CMALIGN_SCORE_FIELDS)
    writer.writerows(read_scores(args.infile))

