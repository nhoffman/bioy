"""
Parse region between primers from fasta file
"""

import logging
import sys

from csv import DictWriter
from itertools import islice, groupby, ifilter

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
    parser.add_argument('-o', '--out-fasta',
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
    ssearch = islice(groupby(parse_ssearch36(args.primer_aligns), lambda hit: hit['q_name']), args.limit)

    aligns = [(q, parse_primer_alignments(h, lprimer = 'lprimer', rprimer = 'rprimer')) for q,h in ssearch]

    aligns = ifilter(lambda (q,h): keep_left(h['l']) and keep_right(h['r']), aligns)

    # a place to put open file handles
    if args.out_rle:
        args.out_rle.writeheader()

    seqs = {f.description:f for f in args.fasta}

    # parse the sequences
    for name, pdict in aligns:
        seq = seqs[name]
        start, stop = pdict['l']['stop'], pdict['r']['start']
        args.out_fasta.write('>{}\n{}\n'.format(seq.description, seq.seq[start:stop]))

        if args.out_rle and args.rle:
            args.out_rle.writerow({'name':seq.name, 'rle':args.rle[seq.name][start:stop]})

