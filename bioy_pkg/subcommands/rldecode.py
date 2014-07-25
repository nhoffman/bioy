"""Run-length decode a fasta file"""

import logging
import sys
import csv

from multiprocessing import Pool

from bioy_pkg.sequtils import homodecode, from_ascii, fastalite
from bioy_pkg import utils
from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('seqs',
                        type=lambda f: fastalite(Opener()(f), readfile=False),
                        help='Input fasta file')
    parser.add_argument('rle',
                        type=Opener(),
                        help='csv file (may be bzip encoded) containing columns "name","rle"')
    parser.add_argument('-o', '--outfile',
                        type=Opener('w'),
                        default=sys.stdout,
                        help='Name of output file')


def seq_and_homodecode(seq_rle):
    seq, rle = seq_rle

    assert len(seq.seq) == len(rle)

    return seq, homodecode(seq.seq, rle)


def action(args):
    # Ignore SIGPIPE, for head support
    utils.exit_on_sigpipe()
    utils.exit_on_sigint()

    rledict = {seqname: rle for seqname, rle in csv.reader(args.rle)}
    seqs = ((s, from_ascii(rledict[s.id])) for s in args.seqs)

    pool = Pool(processes=args.threads)
    seqs = pool.imap(seq_and_homodecode, seqs, chunksize=140)

    for seq, decoded in seqs:
        args.outfile.write('>{}\n{}\n'.format(seq.description, decoded))
