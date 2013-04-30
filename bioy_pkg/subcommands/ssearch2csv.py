"""
Parse ssearch36 -m10 output and print specified contents
"""

import logging
import sys
import pprint
import csv

from itertools import islice, chain, groupby, imap
from operator import itemgetter

from bioy_pkg.sequtils import homodecodealignment, parse_ssearch36, from_ascii
from bioy_pkg.utils import Opener, Csv2Dict

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('alignments',
        default = sys.stdin,
        type = Opener('r'),
        nargs = '?',
        help = 'ssearch -m 10 formatted file')
    parser.add_argument('-o', '--out',
        default = sys.stdout,
        type = Opener('w'),
        help = '(default csv)-formatted output')
    parser.add_argument('-p', '--print-one',
        default = False, action = 'store_true',
        help = 'pretty print first alignment and exit')
    parser.add_argument('-f', '--fieldnames',
        type = lambda f: f.split(','),
        help = 'comma-delimited list of field names to include in output')
    parser.add_argument('--limit',
        type = int,
        metavar = 'N',
        help = 'Print no more than N alignments')
    parser.add_argument('--no-header',
        dest='header',
        action = 'store_false')
    parser.add_argument('-d', '--decode',
        type = Csv2Dict(index = 'name', value = 'rle',
            fieldnames = ['name', 'rle']),
        nargs = '+',
        help = 'Decode alignment')
    parser.add_argument('--min-zscore',
        default = 0,
        type = float,
        metavar = 'X',
        help = 'Exclude alignments with z-score < X')
    parser.add_argument('-a', '--top-alignment',
        default = False, action = 'store_true',
        help = """By default, return all alignments;
                  provide this option to include
                  only the top entry per query.""")

def action(args):
    aligns = islice(parse_ssearch36(args.alignments, False), args.limit)
    aligns = (a for a in aligns if float(a['sw_zscore']) >= args.min_zscore)
    aligns = groupby(aligns, key = itemgetter('q_name'))

    if args.top_alignment:
        aligns = (next(a) for _,a in aligns)
    else:
        aligns = (a for _,i in aligns for a in i) # flatten groupby iters

    if args.decode:
        decoding = {k:v for d in args.decode for k,v in d.items()}
        def decode(aligns):
            aligns['t_seq'], aligns['q_seq'] = homodecodealignment(
                    aligns['t_seq'], from_ascii(decoding[aligns['t_name']]),
                    aligns['q_seq'], from_ascii(decoding[aligns['q_name']]))
            return aligns
        aligns = imap(decode, aligns)

    if args.print_one:
        pprint.pprint(aligns.next())
        sys.exit()

    if args.fieldnames:
        fieldnames = args.fieldnames
    else:
        # peek at first row fieldnames
        top = next(aligns, {})
        fieldnames = top.keys()
        aligns = chain([top], aligns)

    writer = csv.DictWriter(args.out,
            extrasaction = 'ignore',
            fieldnames = fieldnames)

    if args.header:
        writer.writeheader()

    writer.writerows(aligns)

