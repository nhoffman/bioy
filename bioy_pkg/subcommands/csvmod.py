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
Add or rename columns in a csv file.

(rename not yet implemented)
"""

import logging
import sys
import csv
import random
import string
import os
from itertools import imap

from bioy_pkg.utils import Opener, opener, parse_extras

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('infile',
                        default=sys.stdin,
                        type=Opener(),
                        nargs='?',
                        help='input csv file')
    parser.add_argument(
        '-o', '--outfile', default='-', help='output csv file',)
    parser.add_argument('-a', '--add', metavar='FIELD_SPEC',
                        help=("new fields for csv file in form "
                              "'name1:val1,name2:val2,...'"))
    parser.add_argument('-r', '--rename', metavar='FIELD_SPEC',
                        help=("fields to rename in the format "
                              "'from1:to1,from2:to2,...'"))
    parser.add_argument('-i', '--inplace', action='store_true', default=False,
                        help='modify input file in place [%(default)s]')


def tmp(filename):
    dirname, fname = os.path.split(filename)
    return os.path.join(
        dirname,
        ''.join((random.sample(string.letters + string.digits, 15) +
                 ['_' + fname])))


def action(args):
    if args.inplace and args.infile is sys.stdin:
        log.error('Error: cannot use the --inplace option with stdin')
        return

    if args.rename:
        raise NotImplementedError

    reader = csv.DictReader(args.infile)
    fieldnames = reader.fieldnames or []

    new_fields = parse_extras(args.add) if args.add else {}

    if new_fields:
        fieldnames.extend(new_fields.keys())
        reader = imap(lambda row: dict(row, **new_fields), reader)

    if args.inplace:
        outfile = tmp(args.infile.name)
    else:
        outfile = args.outfile

    with opener(outfile, 'w') as fout:
        writer = csv.DictWriter(fout, fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(reader)

    if args.inplace:
        os.rename(fout.name, args.infile.name)
