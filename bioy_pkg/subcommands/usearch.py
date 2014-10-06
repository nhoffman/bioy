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
Run usearch global and produce classify friendly output
"""

import logging
import sys

from itertools import izip, ifilter, imap
from subprocess import Popen, PIPE
from csv import DictWriter

from bioy_pkg.sequtils import BLAST_HEADER, BLAST_FORMAT
from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('fasta',
                        default='-',
                        help='input fasta file')
    parser.add_argument('-o', '--out',
                        type=Opener('w'),
                        default=sys.stdout,
                        help=('tabulated BLAST results with the following headers {}'
                              ).format(BLAST_FORMAT))
    parser.add_argument('-d', '--database',
                        help='a blast database')
    parser.add_argument('--header',
                        action='store_true',
                        help='output header')
    parser.add_argument('--strand',
                        default='plus',
                        choices=['plus', 'minus', 'both'],
                        help="""query strand(s) to search against database/subject.
                      default = %(default)s""")
    parser.add_argument('--id',
                        default=0.9,
                        type=float,
                        help='minimum identity for accepted values default [%(default)s]')
    parser.add_argument('--max',
                        help='maximum number of alignments to keep default = 1')
    parser.add_argument('--usearch', default='usearch6_64',
                        help='name of usearch executable')


def action(args):
    command = [args.usearch]
    command += ['-usearch_global', args.fasta]
    command += ['-threads', str(args.threads)]
    command += ['-id', str(args.id)]
    command += ['-db', args.database]
    command += ['-strand', args.strand]
    command += ['-blast6out', '/dev/stdout']

    if args.max:
        command += ['-maxaccepts', args.max]

    log.debug(' '.join(command))

    usearch = Popen(command, stdout=PIPE, stderr=PIPE)

    lines = imap(lambda l: l.strip().split('\t'), usearch.stdout)

    # usearch has strange commenting at the top it's alignment.
    # we just just want the lines seperated by 12 tabs
    lines = ifilter(lambda l: len(l) == 12, lines)
    lines = imap(lambda l: l[:3] + [l[6], l[7], (int(l[7]) - int(l[6]) + 1)], lines)

    lines = imap(lambda l: izip(BLAST_HEADER, l), lines)
    lines = imap(lambda l: dict(l), lines)

    fieldnames = BLAST_HEADER

    if isinstance(args.coverage, float):
        for l in lines:
            l['coverage'] = (float(l['qend']) - float(l['qstart']) + 1) / float(l['qlen']) * 100
            l['coverage'] = '{0:.2f}'.format(l['coverage'])
        lines = [l for l in lines if float(l['coverage']) >= args.coverage]

        fieldnames += ['coverage']

    out = DictWriter(args.out,
                     fieldnames=BLAST_HEADER,
                     extrasaction='ignore')

    if args.header:
        out.writeheader()

    out.writerows(lines)

    err = usearch.stderr.read().strip()
    if err:
        log.error(err)
