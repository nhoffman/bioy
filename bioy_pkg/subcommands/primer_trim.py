"""
Parse region between primers from fasta file
"""

import logging
import sys

from csv import DictWriter
from itertools import groupby, ifilter, tee
from operator import itemgetter

from bioy_pkg.sequtils import parse_ssearch36, fastalite
from bioy_pkg.utils import Opener, Csv2Dict

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('fasta',
            type = lambda f: fastalite(Opener()(f), readfile = False),
            help = 'input fasta file')
    parser.add_argument('-l', '--left-aligns', type = Opener(),
            help = 'left primer ssearch36 alignment results')
    parser.add_argument('-r', '--right-aligns', type = Opener(),
            help = 'right primer ssearch36 alignment results')
    parser.add_argument('--left-range', metavar = 'START,STOP',
                        help = 'Range of acceptable left primer start positions')
    parser.add_argument('--left-zscore', metavar = 'VALUE', type = float,
                        help = 'Min acceptable left primer z-score')
    parser.add_argument('--right-range', metavar = 'START,STOP',
                        help = 'Range of acceptable right primer start positions')
    parser.add_argument('--right-zscore', metavar = 'VALUE', type = float,
                        help = 'Min acceptable right primer z-score')
    parser.add_argument('--left-expr',
            help = 'python expression defining criteria for keeping left primer')
    parser.add_argument('--right-expr',
            help = 'python expression defining criteria for keeping left primer')
    parser.add_argument('-o', '--fasta-out',
            type = Opener('w'),
            default = sys.stdout,
            help = 'trimmed fasta output file')
    parser.add_argument('--rle',
            type = Csv2Dict('name', 'rle', fieldnames = ['name','rle']),
                        help = 'rle input file (required if --rle-out)')
    parser.add_argument('--rle-out',
            type = lambda f: DictWriter(Opener('w')(f), fieldnames = ['name','rle']),
            help = 'trimmed rle output file')
    parser.add_argument('-i', '--include-primer',
            action = 'store_true', default = False,
            help = 'Include primer in trimmed sequence')

def primer_dict(parsed, side, keep = None, include = False):
    """
    For each alignment between reads (q_seq) and primers (t_seq),
    return a dict of {q_name: position} given whether this is a left
    or right primer (`side`) or whether to include the primer sequence
    (`include`). `position` is a 0-index slice coordinate and is
    defined as follows:

    side   include  -----------------------------------------------------
                        L =======>                       R <=======
    left   False                 ^ (q_al_stop)
    left   True           ^        (q_al_start - 1)
    right  False                                           ^        (q_al_start - 1)
    right  True                                                   ^ (q_al_stop)
    """

    assert side in ('left', 'right')
    assert include in (True, False)

    positions = {
        ('left', False): lambda hit: hit['q_al_stop'],
        ('left', True): lambda hit: hit['q_al_start'] - 1,
        ('right', False): lambda hit: hit['q_al_stop'],
        ('right', True): lambda hit: hit['q_al_start'] - 1
    }

    # 'best' hit is assumed to be first in each group
    hits = (grp.next() for _, grp in groupby(parsed, itemgetter('q_name')))

    if keep:
        hits = ifilter(keep, hits)

    d = {hit['q_name']: int(positions[(side, include)](hit)) for hit in hits}

    msg = '{} sequences passed {} primer criteria'.format(side, len(d))

    if d:
        log.debug(msg)
    else:
        sys.exit(msg)

    return d


def make_fun(expression):
    # I only want eval() to happen once...
    return lambda d: eval(expression)


def make_filter(rangestr, zscore):

    if rangestr and zscore:
        minstart, maxstart = map(int, rangestr.split(','))
        def fun(d):
            start = int(d['q_al_start'])
            return minstart <= start <= maxstart and float(d['sw_zscore']) >= zscore
    elif rangestr:
        minstart, maxstart = map(int, rangestr.split(','))
        def fun(d):
            start = int(d['q_al_start'])
            return minstart <= start <= maxstart
    elif zscore:
        def fun(d):
            return float(d['sw_zscore']) >= zscore
    else:
        raise ValueError('at least one of rangestr and zscore must be provided')

    return fun


def action(args):

    seqs = args.fasta

    left = right = {}

    # right, left are dicts of {name: trim_position}
    if args.right_aligns:
        if args.right_expr:
            keep = make_fun(args.right_expr)
        elif args.right_range or args.right_zscore:
            keep = make_filter(args.right_range, args.right_zscore)
        else:
            keep = None

        right = primer_dict(
            parse_ssearch36(args.right_aligns, numeric = bool(args.right_expr)),
            side='right',
            keep = keep,
            include = args.include_primer)

        seqs = (s for s in seqs if s.id in right)

    if args.left_aligns:
        if args.left_expr:
            keep = make_fun(args.left_expr)
        elif args.left_range or args.left_zscore:
            keep = make_filter(args.left_range, args.left_zscore)
        else:
            keep = None

        left = primer_dict(
            parse_ssearch36(args.left_aligns, numeric = bool(args.left_expr)),
            side='left',
            keep = keep,
            include = args.include_primer)

        seqs = (s for s in seqs if s.id in left)

    seqs, seqs_rle = tee(seqs)

    for s in seqs:
        start, stop = left.get(s.id), right.get(s.id)
        fasta = '>{}\n{}\n'.format(s.description, s.seq[start:stop])
        args.fasta_out.write(fasta)

    if args.rle_out:
        if not args.rle:
            sys.exit('--rle is required')
        args.rle_out.writeheader()

        for s in seqs_rle:
            start, stop = left.get(s.id), right.get(s.id)
            if args.rle_out:
                args.rle_out.writerow({
                    'name': s.id,
                    'rle': args.rle[s.id][start:stop]
                })

