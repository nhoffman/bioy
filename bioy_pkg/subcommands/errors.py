"""
Tally and classify errors given ./ion rlaligns reference and query sequences

Alignments are provided in a csv file with columns ["q_name",
"t_name", "q_align", "t_align"]
"""

import logging
import sys

from collections import Counter
from csv import DictWriter, writer, DictReader
from bioy_pkg import sequtils
from bioy_pkg.sequtils import itemize_errors, error_count
from bioy_pkg.utils import Opener, parse_extras
from os import devnull

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('alignments',
            nargs = '?',
            default = sys.stdin,
            type = Opener('r'),
            help = 'csvfile of errors in format ["q_name", "t_name", "q_align", "t_align"]')
    parser.add_argument('-o', '--out',
            default = sys.stdout,
            type = Opener('w'),
            help = 'csv file tallying each error category for each read')
    parser.add_argument('-m', '--homopolymer-matrix',
            type = Opener('w'),
            default = devnull,
            help = 'csv file containing transition matrix of homopolymer lengths')
    parser.add_argument('-M', '--homopolymer-max',
            default = 6,
            type = int,
            help = 'csv homopolymer length above which counts are binned')
    parser.add_argument('--step',
            action = 'store_true',
            help = 'step through reults (for debugging)')
    parser.add_argument('-f', '--extra-fields',
            help="extra fields for csv file in form 'name1:val1,name2:val2'")
    parser.add_argument('--output-alignment',
            action='store_true',
            help = 'Include the actual alignment in csv output')

def action(args):
    extras = parse_extras(args.extra_fields) if args.extra_fields else {}
    fieldnames = ['t_name', 'q_name', 'length', 'snp',
            'indel', 'homoindel','compound'] + extras.keys()

    if args.output_alignment:
        fieldnames.append('alignment')

    # ssearch query should have been executed with reads as query and
    # reference sequences as target. So each query may have alignments
    # for multiple targets. Therefore, we group by and take only the
    # first (top) hit for each

    # Set up homopolymer matrix vars:
    cnt = Counter()
    gtceil = 'geq{}'.format(args.homopolymer_max)
    ###

    tallies = DictWriter(args.out, fieldnames=fieldnames, extrasaction='ignore')
    tallies.writeheader()

    for a in DictReader(args.alignments):
        e = list(itemize_errors(a['t_seq'], a['q_seq']))

        # instantiate d with zero counts for each error type
        d = {k:0 for k in fieldnames[2:]}
        d['q_name'], d['t_name'] = a['q_name'], a['t_name']
        d.update(error_count(e))
        d.update(extras)

        if args.output_alignment:
            d.update({'alignment':sequtils.show_errors(e)})

        tallies.writerow(d)

        for d in e:
            r, q = d['ref'].strip('=-'), d['query'].strip('=-')
            # count indels or homoindels; exclude compound errors and snps
            if len(set(r+q)) == 1:
                cnt[(len(r) if len(r) <= args.homopolymer_max else gtceil,
                    len(q) if len(q) <= args.homopolymer_max else gtceil)] += 1

        log.debug(a['q_name'])
        log.debug('\n' + sequtils.format_alignment(a['t_seq'], a['q_seq']))
        log.debug(a['q_seq'].replace('-','').replace('=',''))
        log.debug(sequtils.show_errors(e))

        if args.step:
            raw_input()

    # reference counts in rows, query in column
    ii = range(args.homopolymer_max)
    margins = ii + [gtceil]
    w = writer(args.homopolymer_matrix)
    w.writerow(['q{}'.format(i) for i in ii] + [gtceil])
    for i_ref in margins:
        w.writerow(['r{}'.format(i_ref)] + [cnt[i_ref,i_query] for i_query in margins])
