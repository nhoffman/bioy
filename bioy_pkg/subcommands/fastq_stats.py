"""
Describe distributions of sequencing quality scores
"""

import logging
import sys
from csv import DictWriter
from itertools import islice

from numpy import mean
from scipy.stats.mstats import mquantiles
from Bio import SeqIO
from bioy_pkg.utils import Opener, parse_extras

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument(
        'fastq', nargs = '?',
        default = sys.stdin,
        type = Opener('r'),
        help = 'fastq-sanger file with phred scores')
    parser.add_argument(
        '-o', '--out',
        default = sys.stdout,
        type = Opener('w'),
        help = 'csv-format file containing stats for each read.')
    parser.add_argument(
        '-l', '--limit', type = int, metavar = 'N',
        help = 'Limit number of sequences read from input to N')
    parser.add_argument('-e', '--extra-fields',
            help="extra fields for csv file in form 'name1:val1,name2:val2'")
    parser.add_argument('-n', '--no-header', action = 'store_false',
                        default = True, dest = 'show_header')

def action(args):
    extras = parse_extras(args.extra_fields) if args.extra_fields else {}

    quants = [0.025, 0.25, 0.5, 0.75, 0.975]
    fieldnames = ['name', 'length', 'mean']
    fieldnames += ['Q{0:03.0f}'.format(q * 1000) for q in quants]
    fieldnames += extras.keys()

    stats = DictWriter(args.out, fieldnames=fieldnames)
    if args.show_header:
        stats.writeheader()
    for s in islice(SeqIO.parse(args.fastq, 'fastq'), args.limit):
        qual = s.letter_annotations["phred_quality"]
        row = [s.name.replace(':', '_'), len(s), mean(qual)] + \
            list(mquantiles(qual, quants)) + \
            extras.values()
        stats.writerow(dict(zip(fieldnames, row)))
