"""
Calculate the consensus for a multiple aignment
"""

import logging
import sys
import json

from os.path import basename, splitext
from bioy_pkg.sequtils import consensus, fastalite
from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('infile',
            type=Opener(),
            nargs = '?',
            default = sys.stdin,
            help = 'input fasta file (default %(default)s).')
    parser.add_argument('outfile',
            type=Opener('w'),
            default = sys.stdout,
            nargs = '?',
            help='Output fasta file.')
    parser.add_argument('-r', '--rlefile',
            type=Opener(),
            help=('An optional file containing run length '
                  'encoding for infile (.json.bz2)'))
    parser.add_argument('-n', '--seqname',
            help=('Name for the output sequence. '
                  'Default basename(infile).replace(".fasta","") '
                  'if infile is provided, otherwise "consensus"'))
    parser.add_argument('--gaps',
            action = 'store_true',
            help = 'retain gaps in consensus sequence')


def action(args):
    if args.seqname:
        seqname = args.seqname
    else:
        seqname = 'consensus' if args.infile is sys.stdin \
                  else splitext(basename(args.infile.name))[0]

    seqs = list(fastalite(args.infile, 'fasta'))

    if args.rlefile:
        rledict = json.load(args.rlefile)
        rlelist = [rledict[s.id] for s in seqs]
        cons = consensus(seqs, rlelist, degap=not args.gaps)
    else:
        cons = consensus(seqs, degap=not args.gaps)

    args.outfile.write('>{}\n{}\n'.format(seqname, cons))
