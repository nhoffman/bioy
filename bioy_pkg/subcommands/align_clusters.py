"""
Align reads contributing to a denoised cluster.
"""

import logging
import sys
import csv
import re
import subprocess
from itertools import ifilter, groupby
from operator import itemgetter
from os import path
import random

from bioy_pkg.sequtils import SeqLite, fastalite, \
    homodecode, from_ascii, fasta_tempfile
from bioy_pkg.utils import chunker, Opener, Csv2Dict

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('raw_reads',
            type = lambda f: fastalite(Opener()(f)),
            help = """input fasta file containing original
            clustered reads (default stdin).""")
    parser.add_argument('readmap',
            type = Opener('r'),
            help = """output of `bioy denoise --readmap`
            (csv file with columns readname,clustername)""")
    parser.add_argument('-r', '--rlefile',
            type = Csv2Dict('name', 'rle', fieldnames=['name', 'rle']),
            help="""An optional file containing run
            length encoding for infile (.csv.bz2)""")
    parser.add_argument('-d', '--outdir', help='output directory', default='.')
    parser.add_argument('--pattern',
            help = """A regular expression matching cluster names""")
    parser.add_argument('-N', '--sample', type=int, default=100,
                        metavar='N',
                        help='inculde no more than N reads [%(default)s]')


def action(args):
    seqdict = {s.id: s for s in args.raw_reads}

    if args.rlefile:
        def rlemap(seq):
            decoded = homodecode(seq.seq, from_ascii(args.rlefile[seq.id]))
            return SeqLite(seq.id, seq.description, decoded)

    groups = groupby(csv.reader(args.readmap), itemgetter(1))
    for cons, group in groups:
        if args.pattern and not re.search(r'' + args.pattern, cons):
            continue
        reads, _ = zip(*group)
        seqs = [seqdict[name] for name in reads]
        if len(seqs) > args.sample:
            seqs = random.sample(seqs, args.sample)
        if args.rlefile:
            seqs = (rlemap(s) for s in seqs)
            with fasta_tempfile(seqs) as f:
                outfile = path.join(args.outdir, '{}.aln.fasta'.format(cons))
                command = ['muscle', '-quiet', '-seqtype', 'dna',
                           '-in', f, '-out', outfile]
                log.debug(' '.join(command))
                subprocess.check_call(command)
