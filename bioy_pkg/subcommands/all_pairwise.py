"""
Parse ssearch36 -m10 output and print specified contents
"""

import logging
import sys
import pprint
import csv

from itertools import islice, chain, groupby, imap
from operator import itemgetter

from bioy_pkg.sequtils import fasta_tempfile, parse_ssearch36
from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)

def all_pairwise(seqs, cleanup = False):
    """
    Perform all pairwise alignments among sequences in list of SeqRecords `seqs`
    """

    for i in xrange(len(seqs) - 1):
        with fasta_tempfile(seqs[i], 'i%03it_'%i, cleanup) as target:
            with fasta_tempfile(seqs[i+1:], 'i%03iq_'%i, cleanup) as query:
                aligns = run_ssearch36(query, target, cleanup)
                for (d,) in aligns:
                    yield d['t_name'], d['q_name'], d['sw_ident']


def build_parser(parser):
    parser.add_argument('fasta',
        default = sys.stdin,
        type = Opener('r'),
        help = 'sequences in fasta format')
    parser.add_argument('-o', '--out',
        default = sys.stdout,
        type = Opener('w'),
        help = '(default csv)-formatted output')

def action(args):
