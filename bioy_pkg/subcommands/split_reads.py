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

