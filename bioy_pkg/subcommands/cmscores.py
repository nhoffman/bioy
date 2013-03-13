"""
Convert raw cmalign alignment scores to csv format.

Discard the first column ("seq idx") as well as the "elapsed" column
if present.
"""

import logging
import os
import sys
import csv
import argparse
from itertools import imap, groupby, ifilter

from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)

CMALIGN_SCORE_FIELDS = ['name','len','total','struct','prob']

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

    blocks = groupby(infile, lambda line: line.startswith('#'))

    _, comments = blocks.next()
    header, = [line.strip() for line in comments if line.startswith('# seq idx')]

    # Throw out the first column; also throw out the last column if
    # it's called "elapsed" (only present when cmalign is run without
    # mpi).
    if header.endswith('elapsed'):
        rowmap = lambda row: row.split()[1:-1]
    else:
        rowmap = lambda row: row.split()[1:]

    _, data = blocks.next()
    return imap(rowmap, ifilter(lambda row: bool(row.strip()), data))

def build_parser(parser):
    parser.add_argument('infile', type = Opener('r'),
                        help='File containing stdout of call to cmalign')
    parser.add_argument('-o','--outfile', type = Opener('w'),
                        default=sys.stdout,
                        help='Output file in csv format (default is stdout)')

def action(args):

    writer = csv.writer(args.outfile)
    writer.writerow(CMALIGN_SCORE_FIELDS)
    writer.writerows(read_scores(args.infile))
