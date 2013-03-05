"""
Parse barcode, primer, and read from a fastq file
"""

import logging
import os
import csv
from itertools import islice, groupby
import sys

from Bio import SeqIO

from bioy_pkg.sequtils import parse_ssearch36, parse_primer_alignments
from bioy_pkg.utils import Opener, Csv2Dict

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('fasta',
            type = Opener(),
            help = 'input file containing raw reads')
    parser.add_argument('primer_aligns',
            type = Opener(),
            help = 'ssearch36 output (may be bz2-encoded)')
    parser.add_argument('-l','--keep-left',
            help = 'python expression defining criteria for keeping left primer',
            default = "0 <= d.get('start') <= 25 and d.get('sw_zscore') > 50")
    parser.add_argument('-r','--keep-right',
            help = 'python expression defining criteria for keeping left primer',
            default = "200 <= d.get('start') <= 320 and d.get('sw_zscore') > 50")
    parser.add_argument('-d','--outdir',
            help = 'output directory',
            default = '.')
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
            type = Opener('w'),
            default = sys.stdout,
            help = 'trimmed rle output file')

def action(args):
    try:
        os.mkdir(args.outdir)
    except OSError:
        pass

    # I only want eval() to happen once...
    def make_fun(expression):
        return lambda d: eval(expression)

    keep_left = make_fun(args.keep_left)
    keep_right = make_fun(args.keep_right)

    # parse primer alignments
    primerdict = {}
    for q_name, hits in islice(groupby(parse_ssearch36(args.primer_aligns), lambda hit: hit['q_name']), args.limit):
        primerdict[q_name] = parse_primer_alignments(hits, lprimer = 'lprimer', rprimer = 'rprimer')

    # a place to put open file handles
    writer = csv.writer(args.out_rle)
    writer.writerow(['name','rle'])

    # parse the sequences
    for seq in islice(SeqIO.parse(args.fasta, 'fasta'), args.limit):
        pdict = primerdict[seq.name]
        if keep_left(pdict['l']) and keep_right(pdict['r']):
            start, stop = pdict['l']['stop'], pdict['r']['start']
            args.out_fasta.write(
                '>%s\n%s\n' % (seq.name, str(seq.seq)[start:stop]))
            writer.writerow([seq.name, args.rle[seq.name][start:stop]])

