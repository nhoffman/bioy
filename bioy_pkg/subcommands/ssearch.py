"""
Run the ssearch (Smith-Waterman) pairwise aligment tool and output in csv format.

http://computing.bio.cam.ac.uk/local/doc/fasta_guide.pdf

 q_name - query sequence name (ie, from first fasta file)
 t_name - target sequence name (ie, from first fasta file)
 sw_bits - Smith-Waterman bit score
 sw_expect - Smith-Waterman E-value
 sw_frame - alignment frame
 sw_ident - identity score
 sw_overlap - number of overlapping nucleotides
 sw_score - Smith-Waterman score
 sw_sim - similarity score (ie, reduced penalties for ambiguities or similar residues)
 sw_sw_opt - see fasta_guide.pdf section 2.4.2 "Looking at alignments"
 sw_zscore - Smith-Waterman Z-score
 {q,t}_al_display_start - start position of displayed alignment
 {q,t}_al_start - start position  of alignment in sequence
 {q,t}_al_stop - end position  of alignment in sequence
 {q,t}_description - full name from fasta
 {q,t}_seq - nucleotide/AA sequence bounded by {q,t}_al_display_{start,stop}
 {q,t}_sq_len - length of corresponding sequence
 {q,t}_sq_offset - ?
 {q,t}_sq_type - nucleotide or amino acid

Warning: the parse_ssearch36 function does not work with stdin input at this time
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
            help = 'maximum number of alignments to keep [1]')
    parser.add_argument('-g', '--gap-extension-penalty',
            default = '4',
            help = 'gap extension penalty [%(default)s]')
    parser.add_argument('-f', '--gap-open-penalty',
            default = '12',
            help = 'gap open penalty [%(default)s]')
    parser.add_argument('-a', '--full-sequences',
            default = False,
            action = 'store_true',
            help = 'return full sequences in alignment')
    parser.add_argument('--decode',
            type = Csv2Dict(index = 'name', value = 'rle',
                fieldnames = ['name', 'rle']),
            help = 'Decode alignment')
    parser.add_argument('--fieldnames',
            type = lambda f: f.split(','),
            help = 'comma-delimited list of field names to include in output')
    parser.add_argument('-z', '--statistical-calculation',
            default = '1',
            help = """built in statistical calculation
                      of E values for sequences. [%(default)s]""")
    parser.add_argument('--min-zscore',
            default = 0,
            type = float,
            metavar = 'X',
            help = 'Exclude alignments with z-score < X [%(default)s]')

def action(args):
    # setup ssearch command and communicate
    command = ['ssearch36']
    command += ['-m', '10']
    command += ['-3']
    command += ['-n']
    command += ['-z', args.statistical_calculation]
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

    ssearch = Popen(command, stdout = PIPE, stderr = PIPE)

    # parse alignments
    aligns = parse_ssearch36(ssearch.stdout)
    aligns = (a for a in aligns if float(a['sw_zscore']) >= args.min_zscore)
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

    # calculate coverage for each item and repack into generator
    # coverage = |query alignment| / |query length|
    aligns = (dict(d, coverage= str(
             (float(d['q_al_stop']) - float(d['q_al_start']))/float(d['q_sq_len'])
             )) for d in aligns)

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

    error = set(e.strip() for e in ssearch.stderr)
    error = ', '.join(error)

    if ssearch.wait() != 0:
        raise CalledProcessError(ssearch.returncode, error)

    if error:
        log.error(error)
