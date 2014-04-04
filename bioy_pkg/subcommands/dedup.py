"""
Deduplicate sequences by coalescing identical substrings

Use seqmagick convert --deduplicate sequences for fast deduplication
of identical sequences.
"""

import logging
import sys
import csv
from operator import itemgetter

from bioy_pkg.deduplicate import dedup
from bioy_pkg.sequtils import fastalite
from bioy_pkg.utils import Opener, groupbyl

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('sequences',
                        type = Opener(),
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

    if args.split_info:
        primary, secondary = args.primary_group, args.secondary_group
        info_reader = csv.DictReader(args.split_info)
        info = {r['seqname']: r for r in info_reader}

    # group tag sequences if info_file exists
    def group_tag(seq):
        if args.split_info:
            i = info[seq.id]
            group = i[primary] or i[secondary] if secondary else i[primary]
            return dict(group=group, seq=seq)
        else:
            return dict(group=None, seq=seq)

    seqs = (group_tag(s) for s in seqs)

    # group the sequences by tags
    seqs = groupbyl(seqs, key=itemgetter('group'))

    # and now drop the group tags
    seqs = (map(itemgetter('seq'), g) for _,g in seqs)

    # deduplicate each group
    def deduper(seqs):
        """
        Return a list of tuples: (kept_seq, [orig_seqs])
        """

        deduped = dedup([s.seq for s in seqs]).items() # dedup need just the sequences
        deduped = [(seqs[kept], [seqs[o] for o in orig]) for kept,orig in deduped]
        return deduped

    seqs = (deduper(s) for s in seqs)

    # flatten deduped seqs
    seqs = (dp for group in seqs for dp in group)

    # output results
    if args.out_info and args.split_info:
        info_out = csv.DictWriter(args.out_info, fieldnames=info_reader.fieldnames)
        info_out.writeheader()

    if args.out_map:
        map_out = csv.DictWriter(args.out_map, fieldnames = ['kept', 'orig'])

    if args.out_weights:
        weights_out = csv.DictWriter(args.out_weights, fieldnames = ['kept', 'kept', 'weight'])

    for seq,group in seqs:
        args.out.write('>{}\n{}\n'.format(seq.description, seq.seq))

        if args.out_info and args.split_info:
            info_out.writerow(info[seq.id])

        if args.out_map:
            map_out.writerow(dict(kept=seq.id, orig=seq.id))
            for g in group:
                map_out.writerow(dict(kept=seq.id, orig=g.id))

        # add one to include the kept sequence
        if args.out_weights:
            weights_out.writerow(dict(kept=seq.id, weight=len(group) + 1))

