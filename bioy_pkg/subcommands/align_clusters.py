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
    parser.add_argument('--name-suffix',
                        help='string to insert into name before .fasta',
                        default='aln')
    parser.add_argument('--no-align',
                        action='store_false', dest='align', default=True)


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
            log.info(cons)
            reads, _ = zip(*group)
            seqs = [seqdict[name] for name in reads]
            if len(seqs) > args.sample:
                seqs = random.sample(seqs, args.sample)
                if args.rlefile:
                    seqs = (rlemap(s) for s in seqs)
                    outfile = path.join(
                        args.outdir,
                        '{}.{}.fasta'.format(cons, args.name_suffix))

            if args.align:
                with fasta_tempfile(seqs) as f:
                    command = ['muscle', '-quiet', '-seqtype', 'dna',
                               '-in', f, '-out', outfile]
                    log.debug(' '.join(command))
                    subprocess.check_call(command)
            else:
                with open(outfile, 'w') as f:
                    f.write('\n'.join('>{}\n{}'.format(s.id, s.seq) for s in seqs))
