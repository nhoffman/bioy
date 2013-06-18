"""Run-length encode a fasta file"""

import logging
import os
import sys

from itertools import imap, chain
from csv import DictWriter
from multiprocessing import Pool

from bioy_pkg.sequtils import homoencode, to_ascii, fastalite
from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('infiles',
            nargs = '*',
            default = [sys.stdin],
            type = Opener(),
            help = 'Input fasta file')
    parser.add_argument('-o','--outfile',
            type = Opener('w'),
            default = sys.stdout,
            help = 'Name of output file; default: %(default)s')
    parser.add_argument('-r','--rlefile',
            type = lambda f: DictWriter(Opener('w')(f), fieldnames=['name', 'rle']),
            help = """Name of output file for run length encoding; default is to
                      append .csv.bz2 to --outfile basename.""")
    parser.add_argument('--threads',
            default = int(os.environ.get('THREADS_ALLOC') or 1),
            type = int,
            help = """Number of threads (CPUs) to use.
                   Can also specify with environment variable THREADS_ALLOC
                   default = %(default)s""")

def seq_and_homoencode(seq):
    return seq, homoencode(seq.seq)

def action(args):
    pool = Pool(processes = args.threads)

    seqs = imap(fastalite, args.infiles)
    seqs = chain.from_iterable(seqs)
    seqs = pool.imap(seq_and_homoencode, seqs, chunksize = 40)

    for seq, (seqstr, count) in seqs:
        assert len(seqstr) == len(count)
        args.outfile.write('>{}\n{}\n'.format(seq.description, seqstr))
        if args.rlefile:
            args.rlefile.writerow(dict({'name':seq.description, 'rle':to_ascii(count)}))

