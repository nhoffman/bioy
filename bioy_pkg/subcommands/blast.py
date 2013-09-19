"""
Run blastn and produce classify friendly output
"""

import logging
import sys
import os

from itertools import imap, izip
from subprocess import Popen, PIPE
from csv import DictWriter
from cStringIO import StringIO

from bioy_pkg.sequtils import BLAST_FORMAT
from bioy_pkg.utils import Opener

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
    parser.add_argument('--threads',
            default = os.environ.get('THREADS_ALLOC') or '1',
            help = """Number of threads (CPUs) to use in the BLAST search.
                   Can also specify with environment variable THREADS_ALLOC
                   default = %(default)s""")
    parser.add_argument('--id',
            default = '90',
            help = 'minimum identity for accepted values default [%(default)s]')
    parser.add_argument('--max',
            help = 'maximum number of alignments to keep default = (all)')
    parser.add_argument('-n', '--dry-run', action = 'store_true',
                        default = False, help = 'print blast command and exit')

def action(args):
    command = ['blastn']
    command += ['-query', args.fasta]
    command += ['-num_threads', args.threads]
    command += ['-perc_identity', args.id]
    command += ['-outfmt', BLAST_FORMAT]
    command += ['-db', args.database]
    command += ['-strand', args.strand]

    if args.max:
        command += ['-max_target_seqs', args.max]

    if args.dry_run:
        print ' '.join(map(str, command))
        sys.exit(0)

    log.info(command)
    pipe = Popen(command, stdout = PIPE, stderr = PIPE)
    results, errors = pipe.communicate()

    log.error(errors)

    fieldnames = BLAST_FORMAT.split()[1:]

    lines = imap(lambda l: l.strip().split('\t'), StringIO(results))
    lines = imap(lambda l: izip(fieldnames, l), lines)
    lines = imap(lambda l: dict(l), lines)

    out = DictWriter(args.out, fieldnames = fieldnames)

    if args.header:
        out.writeheader()

    out.writerows(lines)
