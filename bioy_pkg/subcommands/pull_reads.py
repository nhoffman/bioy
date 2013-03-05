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
            type = Csv2Dict(fieldnames=['sequence_id','sample_id']),
            help = 'csv(.bz2) file containing sample_id,barcode in the rows.')
    parser.add_argument('-o', '--out',
            type = Opener('w'),
            default = sys.stdout,
            help = 'fasta output file')

def action(args):
    sample_filter = lambda s: args.map_file[s.description] == args.sample_id
    seqs = ifilter(sample_filter, args.fasta)
    args.out.writelines('>{}\n{}\n'.format(s.description, s.seq) for s in seqs)

