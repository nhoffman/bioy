"""
Tally and classify errors given ./ion rlaligns reference and query sequences
"""

import logging
import sys

from collections import Counter
from csv import DictWriter, writer, DictReader
from itertools import imap, tee

from bioy_pkg import sequtils
from bioy_pkg.sequtils import itemize_errors, error_count, show_errors
from bioy_pkg.utils import Opener, parse_extras

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('aligns',
            nargs = '?',
            default = sys.stdin,
            type = Opener('r'),
            help = 'csvfile of ssearch results')
    parser.add_argument('-o', '--out',
            default = sys.stdout,
            type = Opener('w'),
            help = 'csv file tallying each error category for each read')
    parser.add_argument('-m', '--homopolymer-matrix',
            dest = 'matrix',
            type = lambda f: writer(Opener('w')(f)),
            help = 'csv file containing transition matrix of homopolymer lengths')
    parser.add_argument('-M', '--homopolymer-max',
            default = 6,
            type = int,
            help = 'csv homopolymer length above which counts are binned')
    parser.add_argument('--step',
            action = 'store_true',
            help = 'step through reults (for debugging)')
    parser.add_argument('-f', '--extra-fields',
            type = lambda f: parse_extras(f),
            default = {},
            help="extra fields for csv file in form 'name1:val1,name2:val2'")
    parser.add_argument('--output-alignment',
            action='store_true',
            help = 'Include the actual alignment in csv output')

def action(args):
    fieldnames = ['t_name', 'q_name', 'length', 'snp']
    fieldnames += ['indel', 'homoindel','compound']
    fieldnames += args.extra_fields.keys()

    if args.output_alignment:
        fieldnames += ['alignment']

    aligns = DictReader(args.aligns)

    itemizer = lambda a: dict({'errors':itemize_errors(a['t_seq'], a['q_seq'])}, **a)

    aligns = imap(itemizer, aligns)
    aligns, matrix = tee(aligns)

    tallies = DictWriter(args.out,
            fieldnames = fieldnames,
            extrasaction = 'ignore')
    tallies.writeheader()

    # output error counts:
    for a in aligns:
        # instantiate d with zero counts for each error type
        row = {k:0 for k in fieldnames[2:]}
        row['q_name'], row['t_name'] = a['q_name'], a['t_name']
        row.update(error_count(a['errors']))
        row.update(args.extra_fields)

        if args.output_alignment:
            row.update({'alignment':show_errors(a['errors'])})

        tallies.writerow(row)

        log.debug(a['q_name'])
        log.debug('\n' + sequtils.format_alignment(a['t_seq'], a['q_seq']))
        log.debug(a['q_seq'].replace('-','').replace('=',''))
        log.debug(show_errors(a['errors']))

        if args.step:
            raw_input()

    # output homopolymer matrix if specified
    if args.matrix:
        homopolymers = Counter()
        gtceil = 'geq{}'.format(args.homopolymer_max)

        for m in matrix:
            for e in m['errors']:
                r, q = e['ref'].strip('=-'), e['query'].strip('=-')
                # count indels or homoindels; exclude compound errors and snps
                if len(set(r + q)) == 1:
                    ref_count = len(r) if len(r) <= args.homopolymer_max else gtceil
                    query_count = len(q) if len(q) <= args.homopolymer_max else gtceil
                    homopolymers[(ref_count, query_count)] += 1

        # reference counts in rows, query in column
        ii = range(args.homopolymer_max)
        margins = ii + [gtceil]
        args.matrix.writerow(['q{}'.format(i) for i in ii] + [gtceil])
        for i_ref in margins:
            cols = [homopolymers[i_ref, i_query] for i_query in margins]
            args.matrix.writerow(['r{}'.format(i_ref)] + cols)
