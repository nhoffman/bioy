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

"""reverse complement rle and non-rle sequences"""

import logging
import sys
import csv

from bioy_pkg.sequtils import fastalite
from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)

rle_fieldnames = ['name', 'rle']
rev_comp = {'A':'T',
            'T':'A',
            'C':'G',
            'G':'C',
            'M':'K',
            'K':'M',
            'R':'Y',
            'Y':'R',
            'W':'W',
            'S':'S',
            'B':'V',
            'V':'B',
            'D':'H',
            'H':'D',
            'N':'N'}

def build_parser(parser):
    parser.add_argument('infile',
                        type = Opener(),
                        help = 'Input fasta file')
    parser.add_argument('rlefile', nargs = '?', type = Opener(),
                        help = 'csv file (may be bzip encoded) containing columns "name","rle"')
    parser.add_argument('-O', '--out-rle', type = Opener('w'), help = 'reversed rlefile')
    parser.add_argument('-o','--out-fasta',
                        type = Opener('w'),
                        default = sys.stdout,
                        help = 'Name of output file')

def action(args):
    seqs = fastalite(args.infile)

    for s in seqs:
        seq = reversed(s.seq)
        seq = [rev_comp[se] for se in seq]
        seq = ''.join(seq)
        args.out_fasta.write('>{}\n{}\n'.format(s.description, seq))

    if args.rlefile and args.out_rle:
        reader = csv.reader(args.rlefile)
        writer = csv.writer(args.out_rle)

        # try to determine if first row is a header; we'll assume that
        # the first row, second column is a run-length encoding if
        # it's at least half digits.
        name, rle = reader.next()
        if sum(c.isdigit() for c in rle)/float(len(rle)) > 0.5:
            writer.writerow([name, ''.join(reversed(rle))])
        else:
            assert [name, rle] == rle_fieldnames
            writer.writerow([name, rle])

        for name, rle in reader:
            writer.writerow([name, ''.join(reversed(rle))])
