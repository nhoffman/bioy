"""
Add, rename, or remove columns in a csv file.

(rename, remove not yet implemented)
"""

import logging
import sys
import pprint
import csv
import random
import string
import os
from os import path
from itertools import imap

from bioy_pkg.utils import Opener, opener, parse_extras

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('infile',
                        default = sys.stdin,
                        type = Opener('rU'),
                        nargs = '?',
                        help = 'input csv file')
    parser.add_argument('-o', '--outfile', default='-', help = 'output csv file',)
    parser.add_argument('-n', '--new',
                        help="new fields for csv file in form 'name1:val1,name2:val2,...'")
    parser.add_argument('-d', '--delete',
                        help="comma-delimited list of fields to remove")
    parser.add_argument('-r', '--rename',
                        help="fields to rename in the format 'from1:to1,from2:to2,...'")
    parser.add_argument('-i', '--inplace', action='store_true', default=False,
                        help='modify input file in place')


def tmp(filename):
    dirname, fname = path.split(filename)
    return path.join(
        dirname,
        ''.join(random.sample(string.letters + string.digits, 15) + ['_' + fname]))


def action(args):
    if args.inplace and args.infile is sys.stdin:
        sys.exit('Error: cannot use the --inplace option with stdin')

    if args.rename or args.delete:
        raise NotImplementedError

    reader = csv.DictReader(args.infile)
    fieldnames = reader.fieldnames

    new_fields = parse_extras(args.new) if args.new else {}

    if new_fields:
        fieldnames.extend(new_fields.keys())
        reader = imap(lambda row: dict(row, **new_fields), reader)

    with opener(tmp(args.infile.name), 'w') if args.inplace else opener(args.outfile, 'w') as fout:
        writer = csv.DictWriter(fout, fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(reader)

    if args.inplace:
        os.rename(fout.name, args.infile.name)
