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
Calculate all Smith-Waterman pairwise distances among sequences.

TODO: This script needs a lot of work:
      1) The seq_info and tax tables should be
         seperate since they do not exactly match

      2) --split_info should cause --out to output
         [query, target, median]

      3) --all-pairwise should be multithreaded
"""

import logging
import math
import sys
import csv

from numpy import median, around
from operator import itemgetter

from bioy_pkg.sequtils import fastalite, all_pairwise
from bioy_pkg.utils import Opener, groupbyl

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('fasta', default = sys.stdin, type = Opener('r'),
                        nargs = '?', help='sequences in fasta format')
    parser.add_argument('-o', '--out', default = sys.stdout, type = Opener('w'),
                        help='output file (default stdout)')
    parser.add_argument('--matrix-out', type = Opener('w'),
                        help = """median score of pairwise alignments""")
    parser.add_argument('--primary-group', metavar = 'COLUMN_NAME',
                        help = """column in split_info to use for grouping""")
    parser.add_argument('--secondary-group', metavar = 'COLUN_NAME',
                        help = """column in split_info to use for grouping if
                                  primary_group is undefined for a given row""")
    parser.add_argument('--split-info', metavar = 'FILE', type = Opener(),
                        help = """csv file containing column "seqname" plus another
                                  column for grouping sequences prior to deduplication """)
    parser.add_argument('-d', '--distance', action = 'store_true', default = False,
                        help = 'Calculate distance rather than identity.')

def action(args):

    seqs = fastalite(args.fasta)
    pairs = list(all_pairwise(seqs))

    if args.distance:
        pairs = [(q, t, 1 - i) for q, t, i in pairs]

    if args.split_info and args.matrix_out:
        primary, secondary = args.primary_group, args.secondary_group
        split_info = list(csv.DictReader(args.split_info))
        info = {r['seqname']: r for r in split_info if r['seqname']}
        tax = {r['tax_id']:r for r in split_info}

        pairs += map(itemgetter(1,0,2), pairs)

        def group(seqname):
            i = info[seqname]
            return i[primary] or i[secondary] if secondary else i[primary]

        pairs = ((group(left), group(right), score) for left,right,score in pairs)

        # sort and group rows
        pairs = list(groupbyl(pairs, key = itemgetter(0)))

        matrix_out = csv.writer(args.matrix_out)

        # this is the tax_id order we will be using for columns
        tax_ids = map(itemgetter(0), pairs)

        # get the species names to output as first row
        matrix_out.writerow([''] + [tax[t]['tax_name'] for t in tax_ids])

        # iterator through the sorted rows (pairs)
        for row_id, columns in pairs:
            # sort and group columns
            columns = dict(groupbyl(columns, key = itemgetter(1)))

            # get the species name
            row = [tax[row_id]['tax_name']]

            for t in tax_ids:
                # if t not in columns that means there is only
                # sequence representing the group
                # therefore the median destance is 0
                if t not in columns:
                    med = 0
                else:
                    col = columns[t]
                    med = median(map(itemgetter(2), col))
                    # percent and round
                    med = math.ceil(med * 100) / 100

                row.append(med)

            matrix_out.writerow(row)
    else:
        writer = csv.writer(args.out)
        writer.writerow(['query', 'target', 'identity'])
        writer.writerows(pairs)

