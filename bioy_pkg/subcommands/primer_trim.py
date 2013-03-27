"""
Parse region between primers from fasta file
"""

import logging
import sys

from csv import DictWriter
from itertools import islice, groupby, ifilter, imap

from bioy_pkg.sequtils import parse_ssearch36, parse_primer_alignments, fastalite
from bioy_pkg.utils import Opener, Csv2Dict

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('fasta',
            type = lambda f: fastalite(Opener()(f)),
            help = 'input fasta file')
    parser.add_argument('primer_aligns',
            type = Opener(),
            help = 'ssearch36 output (may be bz2-encoded)')
    parser.add_argument('-l','--keep-left',
            help = 'python expression defining criteria for keeping left primer',
            default = "0 <= d.get('start') <= 25 and d.get('sw_zscore') > 50")
    parser.add_argument('-r','--keep-right',
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
            type = Csv2Dict(fieldnames=['name','rle']),
            help = 'rle file')
    parser.add_argument('-O', '--out-rle',
            type = lambda f: DictWriter(Opener('w')(f), fieldnames = ['name','rle']),
            help = 'trimmed rle output file')

def action(args):
    # I only want eval() to happen once...
    def make_fun(expression):
        return lambda d: eval(expression)

    keep_left = make_fun(args.keep_left)
    keep_right = make_fun(args.keep_right)

    # parse primer alignments
    query = lambda hit: hit['q_name']
    aligns = islice(groupby(parse_ssearch36(args.primer_aligns), query), args.limit)

    primer_data = lambda (q,h): (q, parse_primer_alignments(
        h, lprimer = 'lprimer', rprimer = 'rprimer'))
    aligns = imap(primer_data, aligns)

    aligns = ifilter(lambda (q,h): keep_left(h['l']) and keep_right(h['r']), aligns)

    seqs = {f.description:f for f in args.fasta}

    # parse the sequences
    rle_rows = []
    for name, pdict in aligns:
        seq = seqs[name]
        start, stop = pdict['l']['stop'], pdict['r']['start']
        args.out.write('>{}\n{}\n'.format(seq.description, seq.seq[start:stop]))

        if args.rle:
            row = {'name':seq.description, 'rle':args.rle[seq.description][start:stop]}
            rle_rows.append(row)

    if args.out_rle:
        args.out_rle.writeheader()
        args.out_rle.writerows(rle_rows)

