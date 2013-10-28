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
            'V':'B',
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

        header = reader.next()
        assert header == rle_fieldnames

        writer.writerow(header)
        for name, rle in reader:
            writer.writerow([name, ''.join(reversed(rle))])
