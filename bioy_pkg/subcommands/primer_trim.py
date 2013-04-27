"""
Parse region between primers from fasta file
"""

import logging
import sys

from csv import DictWriter
from itertools import groupby, ifilter, imap
from operator import itemgetter

from bioy_pkg.sequtils import parse_ssearch36, fastalite
from bioy_pkg.utils import Opener, Csv2Dict, cast

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('fasta',
            type = lambda f: fastalite(Opener()(f), readfile = False),
            help = 'input fasta file')
    parser.add_argument('-l', '--left',
            type = lambda f: parse_ssearch36(Opener()(f)),
            help = 'left primer ssearch36 alignment results')
    parser.add_argument('-r', '--right',
            type = lambda f: parse_ssearch36(Opener()(f)),
            help = 'right primer ssearch36 alignment results')
    parser.add_argument('--left-expr',
            help = 'python expression defining criteria for keeping left primer',
            default = "0 <= d.get('q_al_start') <= 25 and d.get('sw_zscore') > 50")
    parser.add_argument('--right-expr',
            help = 'python expression defining criteria for keeping left primer',
            default = "200 <= d.get('q_al_start') <= 320 and d.get('sw_zscore') > 50")
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

    # ensure that certain fields are numeric
    def as_numeric(hit):
        return dict(hit, **{k: cast(hit[k]) for k in ['q_al_start', 'q_al_stop', 'sw_zscore']})

    hits = imap(as_numeric, hits)

    if keep:
        hits = ifilter(keep, hits)

    d = {hit['q_name']: positions[(side, include)](hit) for hit in hits}

    msg = '{} sequences passed {} primer criteria'.format(side, len(d))
    if d:
        log.debug(msg)
    else:
        sys.exit(msg)

    return d


def make_fun(expression):
    # I only want eval() to happen once...
    return lambda d: eval(expression)


def action(args):

    seqs = args.fasta

    # right, left are dicts of {name: trim_position}
    if args.right:
        right = primer_dict(
            args.right,
            side='right',
            keep = make_fun(args.right_expr) if args.right_expr else None,
            include = args.include_primer)

    if args.left:
        left = primer_dict(
            args.left,
            side='left',
            keep = make_fun(args.left_expr) if args.left_expr else None,
            include = args.include_primer)

    if args.rle_out:
        if not args.rle:
            sys.exit('--rle is required')
        args.rle_out.writeheader()

    for s in seqs:
        start, stop = 0, len(s)

        # try right first since failure here is more likely
        if args.right:
            try:
                stop = right[s.id]
            except KeyError:
                continue

        if args.left:
            try:
                start = left[s.id]
            except KeyError:
                continue

        fasta = '>{}\n{}\n'.format(s.id, s.seq[start:stop])
        args.fasta_out.write(fasta)

        if args.rle_out:
            args.rle_out.writerow({
                'name': s.id,
                'rle': args.rle[s.id][start:stop]
            })
