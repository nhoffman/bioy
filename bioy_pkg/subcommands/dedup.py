"""
Deduplicate sequences by coalescing identical substrings

Use seqmagick convert --deduplicate sequences for fast deduplication
of identical sequences.
"""

import logging
import sys
import csv
from itertools import groupby, chain, imap
from operator import itemgetter, attrgetter

from bioy_pkg.deduplicate import dedup
from bioy_pkg.sequtils import fastalite
from bioy_pkg.utils import Csv2Dict, opener, Opener

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('sequences',
                        type = lambda f: fastalite(opener(f), readfile = True),
                        help = 'input fasta file')
    parser.add_argument('-i', '--seq-info', metavar='FILE',
                        type = Opener('rU'),
                        help = """csv file containing column "seqname" plus
                      another column for grouping sequences prior to deduplication""")
    parser.add_argument('--primary-group', metavar='COLUMN_NAME',
                        help = 'column in seq_info to use for grouping [default %(default)s]',
                        default = 'species')
    parser.add_argument('--secondary-group', metavar='COLUMN_NAME',
                        help = """column in seq_info to use for grouping
                      if primary_group is undefined for a given row [default %(default)s]""",
                        default = 'tax_id')
    parser.add_argument('-O', '--out-info', metavar='FILE',
                        type = Opener('w'),
                        help = 'deduplicate seq info file')
    parser.add_argument('-o', '--out', metavar='FILE',
                        type = Opener('w'),
                        default = sys.stdout,
                        help = 'deduplicated sequences in fasta format')


def deduper(tups):
    grouper = itemgetter(0)
    tups = sorted(tups, key=grouper)
    for groupname, group in groupby(tups, grouper):
        group = list(group)
        deduped = dedup([s.seq for _, s, _ in group]).keys()
        yield [group[i] for i in deduped]


def action(args):

    if args.seq_info:
        primary, secondary = args.primary_group, args.secondary_group
        info_reader = csv.DictReader(args.seq_info)
        info = {r['seqname']: r for r in info_reader}

        def decorate(seq):
            d = info[seq.id]
            group = d[primary] or d[secondary] if secondary else d[primary]
            return (group, seq, d)

        packed = imap(decorate, args.sequences)
        deduped = chain.from_iterable(deduper(packed))

        if args.out_info:
            info_out = csv.DictWriter(args.out_info, fieldnames=info_reader.fieldnames)
            info_out.writeheader()

        for _, seq, seq_info in deduped:
            args.out.write('>{}\n{}\n'.format(seq.description, seq.seq))
            if args.out_info:
                info_out.writerow(seq_info)
    else:
        seqs = list(args.sequences)
        deduped = dedup([s.seq for s in seqs]).keys()
        args.out.writelines('>{}\n{}\n'.format(
            seqs[i].description, seqs[i].seq) for i in deduped)
