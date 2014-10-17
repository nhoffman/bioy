# This file is part of Bioy
#
#    Bioy is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Bioy is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Bioy.  If not, see <http://www.gnu.org/licenses/>.

"""
Run blastn and produce classify friendly output
"""

import logging
import sys

from csv import DictWriter
from cStringIO import StringIO
from itertools import chain, groupby
from operator import itemgetter
from subprocess import Popen, PIPE

from bioy_pkg.sequtils import BLAST_HEADER, BLAST_FORMAT, fastalite
from bioy_pkg.utils import opener, Opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('fasta',
            default = '-',
            help = 'input fasta file')
    parser.add_argument('-o', '--out',
            type = Opener('w'),
            default = sys.stdout,
            help = 'tabulated BLAST results with the following headers {}'.format(BLAST_HEADER))
    parser.add_argument('-d', '--database', required=True,
            help = 'blast database path for local blasts, or one of "nt" or "nr" if remote')
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
            help = 'minimum identity for accepted values [%(default)s]')
    parser.add_argument('--max',
            help = 'maximum number of alignments to keep default = (all)')
    parser.add_argument('-n', '--dry-run',
            action = 'store_true',
            help = 'print blast command and exit')
    parser.add_argument('--nohits',
            action = 'store_true',
            help = '')
    parser.add_argument('--coverage', type = float,
            help = 'minimum coverage for accepted values [%(default)s]')
    parser.add_argument('--remote', action='store_true',
                        help = 'execute query on remote NCBI server')

def action(args):

    command = ['blastn']
    command += ['-query', args.fasta]
    if args.remote:
        command += ['-remote']
    else:
        command += ['-num_threads', str(args.threads)]
    command += ['-perc_identity', args.id]
    command += ['-outfmt', BLAST_FORMAT]
    command += ['-db', args.database]
    command += ['-strand', args.strand]

    if args.max:
        command += ['-max_target_seqs', args.max]

    log.info(' '.join(command))

    if args.dry_run:
        sys.exit(0)

    pipe = Popen(command, stdout = PIPE, stderr = PIPE)

    results, errors = pipe.communicate()

    if errors:
       log.error(errors)

    # split tab lines
    lines = (r.strip().split('\t') for r in StringIO(results))

    # match with fieldnames
    lines = (zip(BLAST_HEADER, l) for l in lines)

    # make into dict
    lines = [dict(l) for l in lines]

    if isinstance(args.coverage, float):
        for l in lines:
            l['coverage'] = (float(l['qend']) - float(l['qstart']) + 1) \
                    / float(l['qlen']) * 100
            l['coverage'] = '{0:.2f}'.format(l['coverage'])
        lines = [l for l in lines if float(l['coverage']) >= args.coverage]

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

    out = DictWriter(args.out,
                     fieldnames = BLAST_HEADER,
                     extrasaction = 'ignore')

    if args.header:
        out.writeheader()

    out.writerows(lines)
