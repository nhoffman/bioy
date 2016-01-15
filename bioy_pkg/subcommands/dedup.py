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
Fast deduplicate sequences by coalescing identical substrings
"""

import hashlib
import logging
import sys
import csv

from collections import Counter
from itertools import groupby
from operator import itemgetter

from bioy_pkg.sequtils import fastalite
from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('sequences',
                        type = Opener(),
                        default = sys.stdin,
                        help = 'input fasta file')
    parser.add_argument('-i', '--split-info', metavar='FILE',
                        type = Opener('rU'),
                        help = """csv file containing column "seqname" plus
                                  another column for grouping sequences prior to deduplication""")
    parser.add_argument('--primary-group', metavar='COLUMN_NAME',
                        help = 'column in split_info to use for grouping [default %(default)s]',
                        default = 'species')
    parser.add_argument('--secondary-group', metavar='COLUMN_NAME',
                        help = """column in split_info to use for grouping
                                  if primary_group is undefined
                                  for a given row [default %(default)s]""",
                        default = 'tax_id')
    parser.add_argument('--out-map', metavar='FILE',
                        type = Opener('w'),
                        help = 'map file of sequences from (kept_seq_id,orig_seq_i)')
    parser.add_argument('--out-weights', metavar='FILE',
                        help = 'weight file for each kept sequence',
                        type = Opener('w'))
    parser.add_argument('-O', '--out-info', metavar='FILE',
                        type = Opener('w'),
                        help = 'deduplicate seq info file')
    parser.add_argument('-o', '--out', metavar='FILE',
                        type = Opener('w'),
                        default = sys.stdout,
                        help = 'deduplicated sequences in fasta format')

def action(args):
    seqs = fastalite(args.sequences)

    # sort seqs by group information
    if args.split_info:
        primary, secondary = args.primary_group, args.secondary_group
        info_reader = csv.DictReader(args.split_info)
        info = {r['seqname']: r for r in info_reader}

        # group tag sequences if info_file exists
        def group_tag(seq):
            i = info[seq.id]
            group = i[primary] or i[secondary] if secondary else i[primary]
            return dict(group = group, seq = seq)

        seqs = (group_tag(s) for s in seqs)

        # group the sequences by tags
        seqs = sorted(seqs, key = itemgetter('group'))
        seqs = groupby(seqs, key = itemgetter('group'))

        # just need the seqs
        seqs = ((pair['seq'] for pair in group) for _,group in seqs)
    else:
        seqs = (seqs,)

    # set up output files
    if args.out_info and args.split_info:
        info_out = csv.DictWriter(args.out_info, fieldnames=info_reader.fieldnames)
        info_out.writeheader()

    if args.out_map:
        map_out = csv.DictWriter(args.out_map, fieldnames = ['kept', 'orig'])

    if args.out_weights:
        weights_out = csv.DictWriter(args.out_weights, fieldnames = ['kept', 'kept', 'weight'])

    # dedup seqs by groups
    for group in seqs:
        weights = Counter()
        deduped = {}

        for orig in group:
            # checksums are faster to manage
            clean = orig.seq.replace('\n', '').upper()
            checksum = hashlib.sha1(clean).hexdigest()

            if checksum in deduped:
                kept = deduped[checksum]
            else:
                kept = deduped[checksum] = orig
                args.out.write('>{}\n{}\n'.format(kept.description, kept.seq))

                if args.out_info and args.split_info:
                    info_out.writerow(info[kept.id])

            if args.out_weights:
                weights[kept.id] += 1

            if args.out_map:
                map_out.writerow(dict(kept=kept.id, orig=orig.id))

        for kept_id,count in weights.items():
            weights_out.writerow(dict(kept=kept_id, weight=count))


