"""
Run blastn and produce classify friendly output
"""

import logging
import sys
import os

from csv import DictWriter
from cStringIO import StringIO
from itertools import chain, groupby
from operator import itemgetter
from subprocess import Popen, PIPE

from bioy_pkg.sequtils import BLAST_FORMAT, fastalite
from bioy_pkg.utils import opener, Opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('fasta',
            default = '-',
            help = 'input fasta file')
    parser.add_argument('-o', '--out',
            type = Opener('w'),
            default = sys.stdout,
            help = 'tabulated BLAST results with the following headers {}'.format(BLAST_FORMAT))
    parser.add_argument('-d', '--database',
            help = 'a blast database')
    parser.add_argument('--limit',
            type = int,
            help = 'maximum number of query sequences to read from the alignment')
    parser.add_argument('--header',
            action = 'store_true',
            help = 'output header')
    parser.add_argument('--strand',
            default = 'plus',
            choices = ['plus', 'minus', 'both'],
            help = """query strand(s) to search against database/subject.
                      default = %(default)s""")
    parser.add_argument('--id',
            default = '90',
            help = 'minimum identity for accepted values default [%(default)s]')
    parser.add_argument('--max',
            help = 'maximum number of alignments to keep default = (all)')
    parser.add_argument('-n', '--dry-run',
            action = 'store_true',
            help = 'print blast command and exit')
    parser.add_argument('--nohits',
            action = 'store_true',
            help = '')

def action(args):
    command = ['blastn']
    command += ['-query', args.fasta]
    command += ['-num_threads', str(args.threads)]
    command += ['-perc_identity', args.id]
    command += ['-outfmt', BLAST_FORMAT]
    command += ['-db', args.database]
    command += ['-strand', args.strand]

    if args.max:
        command += ['-max_target_seqs', args.max]

    if args.dry_run:
        print ' '.join(map(str, command))
        sys.exit(0)

    log.info(' '.join(command))

    pipe = Popen(command, stdout = PIPE, stderr = PIPE)

    results, errors = pipe.communicate()

    if errors:
       log.error(errors)

    fieldnames = BLAST_FORMAT.split()[1:]

    # split tab lines
    lines = (r.strip().split('\t') for r in StringIO(results))

    # match with fieldnames
    lines = (zip(fieldnames, l) for l in lines)

    # make into dict
    lines = [dict(l) for l in lines]

    if args.nohits:
        # to get nohits first we need to know about the hits
        qids = groupby(lines, key = itemgetter('qseqid'))
        qids = set(q for q,_ in qids)

        # now we can build a list of nohits
        nohits = []
        for q in fastalite(opener(args.fasta)):
            if q.id not in qids:
                nohits.append(q)

        # convert nohits into DictWriter format
        nohits = (dict(qseqid = q.id) for q in nohits)

        # append to lines
        lines = chain(lines, nohits)

    out = DictWriter(args.out, fieldnames = fieldnames)

    if args.header:
        out.writeheader()

    out.writerows(lines)
