"""Run-length decode a fasta file"""

import logging
import os
import sys

from multiprocessing import Pool

from bioy_pkg.sequtils import homodecode, from_ascii, fastalite
from bioy_pkg.utils import Opener, Csv2Dict

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('seqs',
            type = lambda f: fastalite(Opener()(f), readfile = False),
            help = 'Input fasta file')
    parser.add_argument('rle',
            type = Csv2Dict(value = 'rle'),
            help = 'csv file (may be bzip encoded) containing columns "name","rle"')
    parser.add_argument('-o','--outfile',
            type = Opener('w'),
            default = sys.stdout,
            help = 'Name of output file')
    parser.add_argument('--threads',
            default = int(os.environ.get('THREADS_ALLOC') or 1),
            type = int,
            help = """Number of threads (CPUs) to use.
                   Can also specify with environment variable THREADS_ALLOC
                   default = %(default)s""")

def seq_and_homodecode(seq_rle):
    seq, rle = seq_rle

    assert len(seq.seq) == len(rle)

    return seq, homodecode(seq.seq, rle)

def action(args):
    pool = Pool(processes = args.threads)

    seqs = ((s, from_ascii(args.rle[s.id])) for s in args.seqs)
    seqs = pool.imap(seq_and_homodecode, seqs, chunksize = 140)

    for seq,decoded in seqs:
        args.outfile.write('>{}\n{}\n'.format(seq.description, decoded))

