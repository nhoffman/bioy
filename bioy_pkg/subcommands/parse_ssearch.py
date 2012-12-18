"""
Parse ssearch36 -m10 output and print specified contents
"""

import logging
import sys
import pprint
import csv
from itertools import islice

from bioy_pkg.sequtils import homodecodealignment, parse_ssearch36, from_ascii
from bioy_pkg.utils import Opener, csv2dict

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument(
        'alignments', 
        default = sys.stdin,
        type = Opener('r'),
        nargs = '?',
        help = 'ssearch -m 10 formatted file')
    parser.add_argument(
        '-o', '--out', 
        default = sys.stdout,
        type = Opener('w'),
        help = '(default csv)-formatted output')
    parser.add_argument(
        '-p', '--print-one',
        default = False, action = 'store_true',
        help = 'pretty print first alignment and exit')
    parser.add_argument(
        '-f', '--fieldnames',
        default = 'q_name,t_name,sw_ident',
        help = 'comma-delimited list of field names to include in output [%(default)s]')
    parser.add_argument(
            '--all-fieldnames',
            help = 'Output all ssearch fieldnames.  Overrides --fieldnames'
            )
    parser.add_argument(
        '--limit', 
        type = int, 
        metavar = 'N',
        help = 'Print no more than N alignments')
    parser.add_argument(
        '--expr',
        help = """Filter output according to a python expression. 
                  Each alignment is represented by a dict called "a". 
                  Example: --expr="a['sw_zscore'] > 550" """)
    parser.add_argument('--no-header', 
        dest='header', 
        action = 'store_false')
    parser.add_argument('-d', '--decode',
        type = csv2dict('name', 'rle'),
        default = [],
        nargs = '+',
        help = 'Decode alignment')
    parser.add_argument('--fasta', 
            choices = ['t_seq','t_name', 'q_seq', 'q_name', 'all'],
            help = """creates a fasta format of the alignment by 
                      query ('q_seq' or 'q_name'), target ('t_seq' or 't_name') or 'all'""")

def action(args):
    tot = 0
    rledict = {}
    for r in args.decode:
        rledict.update(r)
        tot += len(r)
        log.info('read {} sequences'.format(len(r)))
    # detect name collisions among rle files
    assert len(rledict) == tot

    aligns = islice(parse_ssearch36(args.alignments, False), args.limit)

    if args.print_one:
        pprint.pprint(aligns.next())
        sys.exit()

    if not args.fasta:
        if args.all_fieldnames:
            aligns = list(aligns)
            fieldnames = aligns[0]
        else:
            fieldnames = args.fieldnames.split(',')

        writer = csv.DictWriter(
            args.out, fieldnames=fieldnames,
            extrasaction = 'ignore')

        if args.header:
            writer.writeheader()

    for a in aligns:
        if args.expr:
            if not eval(args.expr):
                continue

        if rledict:
            a['t_seq'], a['q_seq'] = homodecodealignment(
                a['t_seq'], from_ascii(rledict[a['t_name']]),
                a['q_seq'], from_ascii(rledict[a['q_name']]))

        if args.fasta:
            if args.fasta == 'all':
                args.out.write('>{}\n{}\n'.format(a['q_name'], a['q_seq']))
                args.out.write('>{}\n{}\n'.format(a['t_description'], a['t_seq']))
            elif args.fasta in ('q_seq', 'q_name'):
                args.out.write('>{}\n{}\n'.format(a['q_name'], a['q_seq']))
            else:
                args.out.write('>{}\n{}\n'.format(a['t_description'], a['t_seq']))
        else:
            writer.writerow(a)
