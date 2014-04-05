"""
Fast deduplicate sequences by coalescing identical substrings
"""

import hashlib
import logging
import sys
import csv

from collections import OrderedDict
from itertools import chain
from operator import itemgetter

from bioy_pkg.sequtils import fastalite
from bioy_pkg.utils import Opener, groupbyl

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

# deduplicate each group
def deduper(seqs):
    """
    Return a list of tuples: (kept_seq, [orig_seqs])
    """

    deduped = OrderedDict()

    for s in seqs:
        clean = s.seq.replace('\n', '').upper()
        checksum = hashlib.sha1(clean).hexdigest()
        if checksum in deduped:
            deduped[checksum].append(s)
        else:
            deduped[checksum] = [s]

    return [(s.pop(0), s) for s in deduped.values()]

def action(args):
    seqs = fastalite(args.sequences)

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
        seqs = groupbyl(seqs, key = itemgetter('group'))

        # and now drop the group tags
        seqs = (map(itemgetter('seq'), g) for _,g in seqs)

        seqs = chain.from_iterable(deduper(s) for s in seqs)
    else:
        seqs = deduper(seqs)

    # output results
    if args.out_info and args.split_info:
        info_out = csv.DictWriter(args.out_info, fieldnames=info_reader.fieldnames)
        info_out.writeheader()

    if args.out_map:
        map_out = csv.DictWriter(args.out_map, fieldnames = ['kept', 'orig'])

    if args.out_weights:
        weights_out = csv.DictWriter(args.out_weights, fieldnames = ['kept', 'kept', 'weight'])

    for kept,origs in seqs:
        args.out.write('>{}\n{}\n'.format(kept.description, kept.seq))

        if args.out_info and args.split_info:
            info_out.writerow(info[kept.id])

        if args.out_map:
            map_out.writerow(dict(kept=kept.id, orig=kept.id))
            for o in origs:
                map_out.writerow(dict(kept=kept.id, orig=o.id))

        # add one to include the kept sequence
        if args.out_weights:
            weights_out.writerow(dict(kept=kept.id, weight=len(origs) + 1))

