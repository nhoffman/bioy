"""
Return degapped fasta of aligned target region from ssearch results
"""

import logging
import sys

from csv import DictWriter
from itertools import islice, ifilter

from bioy_pkg.sequtils import parse_ssearch36, count_ambiguous
from bioy_pkg.utils import Opener, Csv2Dict

log = logging.getLogger(__name__)

INFO_HEADER = ['seqname', 'tax_id', 'accession', 'description', 'length', 'ambig_count']

def build_parser(parser):
    parser.add_argument('aligns',
        type = lambda f: parse_ssearch36(Opener()(f)),
        help = 'ssearch36 output to pull parsing coordinates')
    parser.add_argument('--info',
        type = Csv2Dict('seqname', fieldnames = INFO_HEADER),
        help = 'seq info file of target seqs')
    parser.add_argument('--min-zscore',
        help = 'minimum z-score',
        type = float)
    parser.add_argument('--rle',
        type = Csv2Dict('name'),
        help = 'run length encoded file (NOT IMPLEMENTED)')
    parser.add_argument('--limit',
        type = int,
        help = 'maximum number of query sequences to read from the alignment')
    parser.add_argument('-O', '--out-info',
        type = lambda f: DictWriter(Opener('w')(f), fieldnames = INFO_HEADER),
        help = 'out info file')
    parser.add_argument('-o', '--out',
        type = Opener('w'),
        default = sys.stdout,
        help = 'output file for fasta reads')

def action(args):
    aligns = ifilter(lambda a: float(a['sw_zscore']) >= args.min_zscore, args.aligns)
    aligns = islice(aligns, args.limit)

    if args.out_info:
        args.out_info.writeheader()

    for a in aligns:
        start = int(a['t_al_start']) - int(a['t_al_display_start'])
        stop = int(a['t_al_stop']) - int(a['t_al_display_start']) + 1
        seq = a['t_seq'].replace('-', '')[start:stop]
        args.out.write('>{}\n{}\n'.format(a['t_description'], seq))

        if args.info and args.out_info:
            seqname = a['t_name']
            info = args.info[seqname]
            info.update({
                'seqname':seqname,
                'length':len(seq),
                'ambig_count':count_ambiguous(seq)
                })
            args.out_info.writerow(info)

