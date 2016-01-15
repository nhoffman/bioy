# This file is part of Bioy
#
#    Bioy is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Bioy is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Bioy.  If not, see <http://www.gnu.org/licenses/>.

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
            args.rlefile.writerow(dict(name = seq.id, rle =  to_ascii(count)))

