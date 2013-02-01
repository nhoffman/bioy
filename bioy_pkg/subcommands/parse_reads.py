"""
Parse barcode, primer, and read from a fastq file
"""

import logging
import os
import csv
from itertools import islice, groupby
import bz2
from os.path import join

from Bio import SeqIO

from bioy_pkg.sequtils import parse_ssearch36, parse_primer_alignments
from bioy_pkg.utils import opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('fasta',
            help = 'input file containing raw reads')
    parser.add_argument('mapfile',
            help = 'csv file containing read name (first column) and specimen (second column); a header row is assumed.')
    parser.add_argument('primer_aligns',
            help = 'ssearch36 output (may be bz2-encoded)')
    parser.add_argument('-l','--keep-left',
            help = 'python expression defining criteria for keeping left primer',
            default = "15 <= d.get('start') <= 25 and d.get('sw_zscore') > 100")
    parser.add_argument('-r','--keep-right',
            help = 'python expression defining criteria for keeping left primer',
            default = "250 <= d.get('start') <= 320 and d.get('sw_zscore') > 100")
    parser.add_argument('-d','--outdir',
            help = 'output directory',
            default = '.')
    parser.add_argument('--limit',
            type = int,
            help = 'maximum number of query sequences to read from the alignment')

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

    # parse barcode mappings
    with opener(args.mapfile) as f:
        bcdict = {row['name']:(row['label'], row['rle']) for row in csv.DictReader(f)}

    # parse primer alignments
    primerdict = {}
    with opener(args.primer_aligns) as f:
        for q_name, hits in islice(groupby(parse_ssearch36(f), lambda hit: hit['q_name']), args.limit):
            primerdict[q_name] = parse_primer_alignments(hits, lprimer = 'lprimer', rprimer = 'rprimer')

    # a place to put open file handles
    fafiles = {}
    rlefiles = {}
    rlewriters = {}
    for sample in set(x[0] for x in bcdict.values()) - {''}:
        fafiles[sample] = open(join(args.outdir,'%s.fasta' % sample),'w')
        rlefile = bz2.BZ2File(join(args.outdir, '%s.csv.bz2' % sample), 'w')
        rlefiles[sample] = rlefile
        rlewriters[sample] = csv.writer(rlefile)
        rlewriters[sample].writerow(['name','rle'])

    # parse the sequences
    with opener(args.fasta) as f:
        for seq in islice(SeqIO.parse(f, 'fasta'), args.limit):
            pdict = primerdict[seq.name]
            if keep_left(pdict['l']) and keep_right(pdict['r']):
                sample, rle = bcdict[seq.name]
                start, stop = pdict['l']['stop'], pdict['r']['start']
                fafiles[sample].write(
                    '>%s\n%s\n' % (seq.name, str(seq.seq)[start:stop]))
                rlewriters[sample].writerow([seq.name, rle[start:stop]])

    # close file objects
    for obj in fafiles.values() + rlefiles.values():
        fname = obj.name
        obj.close()
    #     try:
    #         with opener(fname) as f:
    #             f.next()
    #             f.next()
    #     except StopIteration:
    #         os.remove(fname)
