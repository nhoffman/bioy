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
Parse reads from a fasta file by read to specimen csv map file
"""

import logging

from csv import DictReader
from itertools import groupby
from os import path

from bioy_pkg.sequtils import fastalite
from bioy_pkg.utils import Opener, opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('fasta',
            metavar = 'FILE',
            type = Opener(),
            help = 'input fasta')
    parser.add_argument('specimen_map',
            metavar = 'CSV',
            type = Opener(),
            help = 'columns: readname, specimen')
    parser.add_argument('--outdir',
            metavar = 'DIR',
            help = 'output folder for specimen fasta files, name being specimen.fasta.bz2',
            default = '.')

def action(args):
    fasta = fastalite(args.fasta)

    spec_map = DictReader(args.specimen_map, fieldnames = ['readname', 'specimen'])
    spec_map = {s['readname']:s['specimen'] for s in spec_map}

    def by_specimen(f):
        return spec_map[f.id]

    groups = sorted(fasta, key = by_specimen)
    groups = groupby(groups, key = by_specimen)

    for spec, fasta in groups:
        fasta = ('>{}\n{}'.format(f.description, f.seq) for f in fasta)
        fasta = '\n'.join(fasta)

        filename = path.join(args.outdir, '{}.fasta.bz2'.format(spec))

        with opener(filename, 'w') as out:
            out.write(fasta)

