"""
Run the fasta pairwise aligment tool and output in csv format.

http://computing.bio.cam.ac.uk/local/doc/fasta_guide.pdf
"""

import logging
import sys

from itertools import chain, groupby, imap, islice
from operator import itemgetter
from cStringIO import StringIO
from subprocess import Popen, PIPE
from csv import DictWriter

from bioy_pkg.sequtils import fastalite, parse_ssearch36, homodecodealignment, from_ascii
from bioy_pkg.utils import Opener, Csv2Dict

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('query',
            default = sys.stdin,
            type = Opener(),
            help = 'input fasta query file')
    parser.add_argument('-l', '--library',
            required = True,
            help = 'input fasta library file to search against')
    parser.add_argument('-o', '--out',
            type = Opener('w'),
            default = sys.stdout,
            help = 'tabulated ssearch results')
    parser.add_argument('--limit',
            type = int,
            help = 'maximum number of query sequences to read from the alignment')
    parser.add_argument('--no-header',
            dest = 'header',
            action = 'store_false',
            default = True,
            help = 'no header')
    parser.add_argument('--all-alignments',
            action = 'store_true',
            help = 'maximum number of alignments to keep default = 1')
    parser.add_argument('-g', '--gap-extension-penalty',
            default = '4',
            help = 'gap extension penalty default = %(default)s')
    parser.add_argument('-f', '--gap-open-penalty',
            default = '12',
            help = 'gap open penalty default = %(default)s')
    parser.add_argument('-a', '--full-sequences',
            default = False,
            action = 'store_true',
            help = 'return full sequences in alignment')
    parser.add_argument('-O', '--out-raw',
            type = Opener('w'),
            help = 'return raw ssearch output')
    parser.add_argument('--decode',
            type = Csv2Dict(index = 'name', value = 'rle',
                fieldnames = ['name', 'rle']),
            help = 'Decode alignment')
    parser.add_argument('--fieldnames',
            type = lambda f: f.split(','),
            help = 'comma-delimited list of field names to include in output')
    parser.add_argument('--min-zscore',
            default = 0,
            type = float,
            metavar = 'X',
            help = 'Exclude alignments with z-score < X')

def action(args):
    # setup ssearch command and communicate
    command = ['fasta36']
    command += ['-m', '10']
    command += ['-3']
    command += ['-n']
    command += ['-g', args.gap_extension_penalty]
    command += ['-f', args.gap_open_penalty]
    command += ['-T', args.threads]

    if args.full_sequences:
        command += ['-a']

    if not args.all_alignments:
        command += ['-b', '1']
        command += ['-d', '1']

    command += ['@', args.library]

    pipe = Popen(command, stdout = PIPE, stderr = PIPE, stdin = PIPE)

    fasta = fastalite(args.query, readfile = False)
    fasta = islice(fasta, args.limit)
    fasta = ('>{}\n{}\n'.format(f.description, f.seq) for f in fasta)
    fasta = ''.join(fasta)

    results, errors = pipe.communicate(fasta)

    log.error(errors)

    if args.out_raw:
        args.out_raw.write(results)

    # parse alignments
    aligns = StringIO(results)
    aligns = parse_ssearch36(aligns)
    aligns = (a for a in aligns if float(a['fa_zscore']) >= args.min_zscore)
    aligns = groupby(aligns, key = itemgetter('q_name'))
    aligns = (a for _,i in aligns for a in i) # flatten groupby iters

    # decode if appropriate
    if args.decode:
        decoding = {k:v for d in args.decode for k,v in d.items()}
        def decode(aligns):
            aligns['t_seq'], aligns['q_seq'] = homodecodealignment(
                    aligns['t_seq'], from_ascii(decoding[aligns['t_name']]),
                    aligns['q_seq'], from_ascii(decoding[aligns['q_name']]))
            return aligns
        aligns = imap(decode, aligns)

    # write results
    if args.fieldnames:
        fieldnames = args.fieldnames
    else:
        # peek at first row fieldnames
        top = next(aligns, {})
        fieldnames = top.keys()
        aligns = chain([top], aligns)

    writer = DictWriter(args.out,
            extrasaction = 'ignore',
            fieldnames = fieldnames)

    if args.header:
        writer.writeheader()

    writer.writerows(aligns)

