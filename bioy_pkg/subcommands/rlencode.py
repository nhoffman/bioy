"""Run-length encode a fasta file"""

import logging
import sys

from itertools import imap, chain
from csv import DictWriter
from multiprocessing import Pool

from bioy_pkg.sequtils import homoencode, to_ascii, fastalite
from bioy_pkg import utils
from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('infiles',
                        nargs='*',
                        default=[sys.stdin],
                        type=Opener(),
                        help='Input fasta file')
    parser.add_argument('-o', '--outfile',
                        type=Opener('w'),
                        default=sys.stdout,
                        help='Name of output file; default: %(default)s')
    parser.add_argument('-r', '--rlefile',
                        type=lambda f: DictWriter(
                            Opener('w')(f), fieldnames=['name', 'rle']),
                        help="""Name of output file for run length encoding; default is to
                      append .csv.bz2 to --outfile basename.""")


def seq_and_homoencode(seq):
    return seq, homoencode(seq.seq)


def action(args):
    # Ignore SIGPIPE, for head support
    utils.exit_on_sigpipe()
    utils.exit_on_sigint()

    pool = Pool(processes=args.threads)

    seqs = imap(fastalite, args.infiles)
    seqs = chain.from_iterable(seqs)
    seqs = pool.imap(seq_and_homoencode, seqs, chunksize=40)

    for seq, (seqstr, count) in seqs:
        assert len(seqstr) == len(count)

        args.outfile.write('>{}\n{}\n'.format(seq.description, seqstr))

        if args.rlefile:
            args.rlefile.writerow(dict(name=seq.id, rle=to_ascii(count)))
