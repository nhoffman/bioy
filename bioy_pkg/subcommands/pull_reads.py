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
Parse barcode, primer, and read from a fastq file
"""

import logging
import sys

from itertools import ifilter

from bioy_pkg.sequtils import fastalite
from bioy_pkg.utils import opener, Opener, Csv2Dict

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('fasta',
            type = lambda f: fastalite(opener(f)),
            help = 'input file containing raw reads')
    parser.add_argument('--sample-id',
            help = 'sample id to pull reads for')
    parser.add_argument('--map-file',
            type = Csv2Dict(value = 'sample_id', fieldnames=['sequence_id','sample_id']),
            help = 'csv(.bz2) file containing sequence_id,sample_id in the rows.')
    parser.add_argument('-o', '--out',
            type = Opener('w'),
            default = sys.stdout,
            help = 'fasta output file')

def action(args):
    sample_filter = lambda s: args.map_file[s.description] == args.sample_id
    seqs = ifilter(sample_filter, args.fasta)
    args.out.writelines('>{}\n{}\n'.format(s.description, s.seq) for s in seqs)

