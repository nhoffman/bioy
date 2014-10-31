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
import csv

from itertools import ifilter
from subprocess import Popen, PIPE

from bioy_pkg.sequtils import BLAST_HEADER, USEARCH_BLAST6OUT_HEADERS, BLAST_FORMAT
from bioy_pkg.utils import Opener, named_tempfile

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
    parser.add_argument('--min-coverage', type=float,
                        help='minimum percent coverage for each alignment [%(default)s]')
    parser.add_argument('--max',
                        help='maximum number of alignments to keep default = 1')
    parser.add_argument('--usearch', default='usearch6_64',
                        help='name of usearch executable')


def parse_usearch(lines):
    """Return an iterable of dicts from output of 'usearch -blast6out'. See
    http://drive5.com/usearch/manual/blast6out.html for output format.

    Coverage is calculated relative to the length of the query sequence.
    """

    fieldnames = USEARCH_BLAST6OUT_HEADERS[:]
    fieldnames[fieldnames.index('length')] = 'qlen'
    reader = csv.DictReader(lines, fieldnames=fieldnames, delimiter='\t')

    for d in reader:
        d['coverage'] = 100.0 * (float(d['qend']) - float(d['qstart']) + 1) / float(d['qlen'])

        yield d


def test_parse_usearch():

    lines = """Actinomyces|6|IBRIB9O01DNL9H	S003710619	99.1	442	4	0	1	442	1	1501	*	*
Actinomyces|2|IBRIB9O01B0977	S002449772	98.9	445	3	0	1	443	1	1394	*	*
Actinomyces|3|IBRIB9O01AV846	S002952986	99.1	447	2	0	1	447	1	1377	*	*""".splitlines()

    result = list(parse_usearch(lines))
    assert len(result) == 3
    for k in BLAST_HEADER:
        assert k in result[0]


def action(args):

    with named_tempfile('rw') as tfile:
        command = [args.usearch,
                   '-usearch_global', args.fasta,
                   '-threads', str(args.threads),
                   '-id', str(args.id),
                   '-db', args.database,
                   '-strand', args.strand,
                   '-blast6out', tfile.name]

        if args.max:
            command += ['-maxaccepts', args.max]

        log.info(' '.join(command))

        usearch_proc = Popen(command, stderr=PIPE, stdout=PIPE)

        errmsg = usearch_proc.stderr.read()

        usearch_proc.communicate()
        if usearch_proc.returncode != 0:
            log.error(errmsg)
            return usearch_proc.returncode

        tfile.flush()
        tfile.seek(0)
        results = parse_usearch(tfile)

        if args.min_coverage:
            results = ifilter(lambda d: d['coverage'] >= args.min_coverage, results)

        writer = csv.DictWriter(args.out, fieldnames=BLAST_HEADER, extrasaction='ignore')

        if args.header:
            writer.writeheader()

        writer.writerows(results)
