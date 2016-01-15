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

"""Run-length decode a fasta file"""

import logging
import os
import sys
import csv

from multiprocessing import Pool

from bioy_pkg.sequtils import homodecode, from_ascii, fastalite
from bioy_pkg.utils import Opener, Csv2Dict

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('seqs',
            type = lambda f: fastalite(Opener()(f), readfile = False),
            help = 'Input fasta file')
    parser.add_argument('rle',
            type = Opener(),
            help = 'csv file (may be bzip encoded) containing columns "name","rle"')
    parser.add_argument('-o','--outfile',
            type = Opener('w'),
            default = sys.stdout,
            help = 'Name of output file')

def seq_and_homodecode(seq_rle):
    seq, rle = seq_rle

    assert len(seq.seq) == len(rle)

    return seq, homodecode(seq.seq, rle)

def action(args):

    rledict = {seqname: rle for seqname, rle in csv.reader(args.rle)}
    seqs = ((s, from_ascii(rledict[s.id])) for s in args.seqs)

    pool = Pool(processes = args.threads)
    seqs = pool.imap(seq_and_homodecode, seqs, chunksize = 140)

    for seq,decoded in seqs:
        args.outfile.write('>{}\n{}\n'.format(seq.description, decoded))

