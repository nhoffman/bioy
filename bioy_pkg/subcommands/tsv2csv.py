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

import csv
import logging
import re
import sys

from functools import partial

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
            default = sys.stdout,
            help = 'csv file')
    parser.add_argument('--split-column',
            help = 'column:delimiter:*newcolumns')

def parse_column(column):
    colons = re.compile(r"""((?:[^:"']|"[^"]*"|'[^']*')+)""")
    return colons.split(column)[1::2]

def action(args):
    inn = csv.DictReader(args.csv, dialect=csv.excel_tab)
    fieldnames = inn.fieldnames

    if args.split_column:
        parsed = parse_column(args.split_column)

        if len(parsed) < 3:
            sys.exit('please enter at least three columns "column:delimiter:*newcolumns"')

        (column, delim), newnames = parsed[:2], parsed[2:]
        index = fieldnames.index(column) + 1
        fieldnames = fieldnames[:index] + newnames + fieldnames[index:]

        def make_new_columns(line):
            newcolumns = line[column].split(delim)

            if len(newcolumns) == len(newnames):
                newcolumns = [c.strip() for c in newcolumns]
                newcolumns = zip(newnames, newcolumns)
                return dict(newcolumns, **line)
            else:
                return line

        inn = (make_new_columns(l) for l in inn)

    out = csv.DictWriter(args.out, dialect=csv.excel, fieldnames=fieldnames)
    out.writeheader()

    for l in inn:
        out.writerow(l)

