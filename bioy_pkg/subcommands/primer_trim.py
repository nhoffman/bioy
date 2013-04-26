"""
Parse region between primers from fasta file
"""

import logging
import sys

from csv import DictWriter
from itertools import groupby, ifilter, imap, tee
from operator import itemgetter

from bioy_pkg.sequtils import parse_ssearch36, fastalite
from bioy_pkg.utils import Opener, Csv2Dict

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('fasta',
            type = lambda f: fastalite(Opener()(f), readfile = False),
            help = 'input fasta file')
    parser.add_argument('-l', '--left',
            type = lambda f: parse_ssearch36(Opener()(f)),
            default = [],
            help = 'left primer ssearch36 alignment results')
    parser.add_argument('-r', '--right',
            type = lambda f: parse_ssearch36(Opener()(f)),
            default = [],
            help = 'right primer ssearch36 alignment results')
    parser.add_argument('--left-expr',
            help = 'python expression defining criteria for keeping left primer',
            default = "0 <= d.get('start') <= 25 and d.get('sw_zscore') > 50")
    parser.add_argument('--right-expr',
            help = 'python expression defining criteria for keeping left primer',
            default = "200 <= d.get('start') <= 320 and d.get('sw_zscore') > 50")
    parser.add_argument('--limit',
            type = int,
            help = 'maximum number of query sequences to read from the alignment')
    parser.add_argument('-o', '--out',
            type = Opener('w'),
            default = sys.stdout,
            help = 'trimmed fasta output file')
    parser.add_argument('--rle',
            type = Csv2Dict('name', 'rle', fieldnames = ['name','rle']),
            help = 'rle file')
    parser.add_argument('-O', '--out-rle',
            type = lambda f: DictWriter(Opener('w')(f), fieldnames = ['name','rle']),
            help = 'trimmed rle output file')

def simple_data(hit):
    d = {}

    d['start']= int(hit['q_al_start']) - 1
    d['stop'] = int(hit['q_al_stop'])
    d['sw_zscore'] = float(hit['sw_zscore'])
    d['q_name'] = hit['q_name']

    return d

def action(args):
    # I only want eval() to happen once...
    def make_fun(expression):
        return lambda d: eval(expression)

    keep_left = make_fun(args.left_expr)
    keep_right = make_fun(args.right_expr)

    left = right = {}

    if args.left:
        left = (simple_data(l) for l in args.left)
        left = groupby(left, itemgetter('q_name'))
        left = imap(lambda (q,h):
                (q, sorted(h, key = itemgetter('sw_zscore'), reverse = True)[0]), left)
        left = ifilter(lambda (q,h): keep_left(h), left)
        left = dict(left)

    if args.right:
        right = ifilter(lambda r: r['q_name'] in left, args.right) if left else args.right
        right = (simple_data(r) for r in right)
        right = groupby(right, itemgetter('q_name'))
        right = imap(lambda (q,h):
                (q, sorted(h, key = itemgetter('sw_zscore'), reverse = True)[0]), right)
        right = ifilter(lambda (q,h): keep_right(h), right)
        right = dict(right)

    seqs = (s for s in args.fasta if s.description in (right or left))
    seqs, rle = tee(seqs)

    for s in seqs:
        start = left.get('stop', None)
        stop = right.get('start', None)
        fasta = '>{}\n{}\n'.format(s.description, s.seq[start:stop])
        args.out.write(fasta)

    # parse the rle's
    if args.out_rle:
        args.out_rle.writeheader()
        for s in rle:
            start = left.get('stop', None)
            stop = right.get('start', None)
            row = {'name':s.description,
                   'rle':args.rle[s.description][start:stop]}
            args.out_rle.writerow(row)

