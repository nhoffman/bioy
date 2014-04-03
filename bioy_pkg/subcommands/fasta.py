"""
Run the fasta pairwise aligment tool and output in csv format.

http://computing.bio.cam.ac.uk/local/doc/fasta_guide.pdf
"""

import logging
import sys

from itertools import chain, groupby, imap
from operator import itemgetter
from subprocess import Popen, PIPE, CalledProcessError
from csv import DictWriter

from bioy_pkg.sequtils import parse_ssearch36, homodecodealignment, from_ascii
from bioy_pkg.utils import Opener, Csv2Dict

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('query',
            help = 'input fasta query file')
    parser.add_argument('library',
            help = 'input fasta library file to search against')
    parser.add_argument('-o', '--out',
            type = Opener('w'),
            default = sys.stdout,
            help = 'tabulated ssearch results')
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
    command += ['-T', str(args.threads)]

    if args.full_sequences:
        command += ['-a']

    if not args.all_alignments:
        command += ['-b', '1']
        command += ['-d', '1']

    command += [args.query, args.library]

    log.info(' '.join(command))

    pipe = Popen(command, stdout = PIPE, stderr = PIPE)

    # parse alignments
    aligns = parse_ssearch36(pipe.stdout)
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
        if top:
            aligns = chain([top], aligns)

    if fieldnames:
        writer = DictWriter(args.out,
                extrasaction = 'ignore',
                fieldnames = fieldnames)

        if args.header:
            writer.writeheader()

        for a in aligns:
            writer.writerow(a)

    error = set(e.strip() for e in pipe.stderr)
    error = ', '.join(error)

    if pipe.wait() != 0:
        raise CalledProcessError(pipe.returncode, error)
    if error:
        log.error(error)

