"""
Parse region between primers from fasta file
"""

import logging
import sys

from csv import DictWriter
from itertools import islice, groupby, ifilter, imap, tee, chain
from operator import itemgetter

from bioy_pkg.sequtils import parse_ssearch36, parse_primer_alignments, fastalite
from bioy_pkg.utils import Opener, Csv2Dict

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('fasta',
            type = lambda f: fastalite(Opener()(f), readfile = False),
            help = 'input fasta file')
    parser.add_argument('primer_aligns',
            nargs = '+',
            type = Opener(),
            help = 'ssearch36 output (may be bz2-encoded)')
    parser.add_argument('-l','--keep-left',
            help = 'python expression defining criteria for keeping left primer',
            default = "0 <= d.get('start') <= 25 and d.get('sw_zscore') > 50")
    parser.add_argument('-r','--keep-right',
            help = 'python expression defining criteria for keeping left primer',
            default = "200 <= d.get('start') <= 320 and d.get('sw_zscore') > 50")
    parser.add_argument('--left-primer-name',
            default = 'lprimer',
            help = 'name of left primer. default = %(default)s')
    parser.add_argument('--right-primer-name',
            default = 'rprimer',
            help = 'name of right primer. default = %(default)s')
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

def action(args):
    # I only want eval() to happen once...
    def make_fun(expression):
        return lambda d: eval(expression)

    keep_left = make_fun(args.keep_left)
    keep_right = make_fun(args.keep_right)

    # apply custom primer function to each alignment
    primer_data = lambda (q,h): (q, parse_primer_alignments(
        h, lprimer = args.left_primer_name, rprimer = args.right_primer_name))

    # combine ssearch aligns
    aligns = (parse_ssearch36(s) for s in args.primer_aligns)
    aligns = chain(*aligns)
    aligns = sorted(aligns, key = itemgetter('q_name'))
    aligns = groupby(aligns, itemgetter('q_name'))
    aligns = islice(aligns, args.limit)
    aligns = imap(primer_data, aligns) # parse primer information
    aligns = ifilter(lambda (q,h): keep_left(h['l']), aligns)
    aligns = ifilter(lambda (q,h): keep_right(h['r']), aligns)

    aligns, aligns_rle = tee(aligns)

    seqs = {f.id:f for f in args.fasta}

    # parse the sequences
    for name, pdict in aligns:
        seq = seqs[name]
        start, stop = pdict['l']['stop'], pdict['r']['start']
        fasta = '>{}\n{}\n'.format(seq.description, seq.seq[start:stop])
        args.out.write(fasta)

    # parse the rle's
    if args.out_rle:
        args.out_rle.writeheader()
        for name, pdict in aligns_rle:
            seq = seqs[name]
            row = {'name':seq.description,
                   'rle':args.rle[seq.description][start:stop]}
            args.out_rle.writerow(row)

