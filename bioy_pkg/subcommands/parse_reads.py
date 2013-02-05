"""
Parse barcode, primer, and read from a fastq file
"""

import logging
import csv
import sys

from itertools import islice, groupby
from os.path import join, basename, splitext, dirname

from bioy_pkg.sequtils import parse_ssearch36, parse_primer_alignments, fastalite
from bioy_pkg.utils import opener, Opener, Csv2Dict

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('fasta',
            type = Opener(),
            help = 'input file containing raw reads')
    parser.add_argument('--mapfile',
            type = Csv2Dict('name'),
            help = 'csv file containing read name (first column) and specimen (second column); a header row is assumed. Overwrites --rle-file, --out')
    parser.add_argument('--primer-aligns',
            type = Opener(),
            help = 'ssearch36 output (may be bz2-encoded)')
    parser.add_argument('-l','--keep-left',
            help = 'python expression defining criteria for keeping left primer',
            default = "d.get('start') <= 100 and d.get('sw_zscore') > 100")
    parser.add_argument('-r','--keep-right',
            help = 'python expression defining criteria for keeping left primer',
            default = "200 <= d.get('start') and d.get('sw_zscore') > 100")
    parser.add_argument('--rle',
            type = Csv2Dict('name'),
            help = 'run length encoded file')
    parser.add_argument('--limit',
            type = int,
            help = 'maximum number of query sequences to read from the alignment')
    parser.add_argument('-o', '--out',
            help = """output file for fasta reads.  If not mapfile is specified the
                      basename of this file will serve as the basename for the rle file as well.""")

def action(args):
    # I only want eval() to happen once...
    def make_fun(expression):
        return lambda d: eval(expression)

    keep_left = make_fun(args.keep_left)
    keep_right = make_fun(args.keep_right)

    # parse primer alignments
    primerdict = {}
    for q_name, hits in islice(groupby(parse_ssearch36(args.primer_aligns),
        lambda hit: hit['q_name']), args.limit):
        primerdict[q_name] = parse_primer_alignments(hits, lprimer = 'lprimer', rprimer = 'rprimer')

    # a place to put open file handles
    fafiles = {}
    rlewriters = {}
    if args.mapfile:
        for sample in set(x['label'] for x in args.mapfile.values()) - {''}:
            fafiles[sample] = opener(join(args.outdir,'{}.fasta'.format(sample)), 'w')
            rlefile = opener(join(args.outdir, '{}.csv.bz2'.format(sample)), 'w')
            rlewriters[sample] = csv.writer(rlefile)
            rlewriters[sample].writerow(['name','rle'])
    else:
        fa_base,_ = splitext(args.out)
        fa_name = basename(fa_base)
        dir_name = dirname(fa_base)
        fafiles[fa_name] = opener(join(dir_name, '{}.fasta'.format(fa_name)), 'w')
        rlewriters[fa_name] = csv.writer(opener(join(dir_name, '{}.csv.bz2'.format(fa_name)), 'w'))
        rlewriters[fa_name].writerow(['name','rle'])

    # parse the sequences
    for seq in islice(fastalite(args.fasta, readfile = False), args.limit):
        pdict = primerdict[seq.description]
        start = pdict['l']['stop'] if keep_left(pdict['l']) else 0
        stop = pdict['r']['start'] if keep_right(pdict['r']) else len(seq.seq)

        if args.mapfile:
            sample, rle = args.mapfile[seq.description]
        else:
            sample, rle = fa_name, args.rle[seq.description]

        fafiles[sample].write('>{}\n{}\n'.format(seq.description, str(seq.seq)[start:stop]))
        rlewriters[sample].writerow([seq.description, rle[start:stop]])

