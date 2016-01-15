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
import os
import sys
import csv

from subprocess import Popen, PIPE, CalledProcessError

from bioy_pkg.utils import Opener, named_tempfile
from bioy_pkg.sequtils import USEARCH_HEADER

log = logging.getLogger(__name__)

toSsearch = {'qseqid': 'q_name',
             'sseqid': 't_name',
             'pident': 'sw_ident',
             'length': 'q_sq_len',
             'mismatch': 'mismatch',
             'qcovs': 'qcovs',
             'gapopen': 'gaopen',
             'qstart': 'q_al_start',
             'qend': 'q_al_stop',
             'sstart': 't_al_start',
             'send': 't_al_stop',
             'evalue': 'sw_expect',
             'bitscore': 'sw_bits'}


def build_parser(parser):
    parser.add_argument('query',
                        help='input fasta query file')
    parser.add_argument('library',
                        help='input fasta library file to search against')
    parser.add_argument('-o', '--out',
                        type=Opener('w'),
                        default=sys.stdout,
                        help='tabulated ssearch results')
    parser.add_argument('--no-header',
                        dest='header',
                        action='store_false',
                        help='no header')
    parser.add_argument('--strand',
                        default='plus',
                        choices=['plus', 'minus', 'both'],
                        help=('query strand(s) to search against '
                              'database/subject.default = %(default)s'))
    parser.add_argument('--min-identity',
                        default=0.9,
                        type=float,
                        help=('minimum identity for accepted '
                              'values default [%(default)s]'))
    parser.add_argument('--min-coverage', type=float,
                        help='minimum percent coverage for each alignment')
    parser.add_argument('--max',
                        help=('maximum number of alignments '
                              'to keep [%(default)s]'))
    parser.add_argument('--usearch', default='usearch6',
                        help='name of usearch executable')
    parser.add_argument('--fieldnames',
                        help='specify fieldnames to use')


def coverage(d):
    """
    Coverage is calculated relative to the length of the query sequence.
    """

    d['qcovs'] = float(d['qend']) - float(d['qstart']) + 1
    d['qcovs'] /= float(d['length'])

    return d


def sw_ident_to_decimal(result):
    result['sw_ident'] = float(result['sw_ident']) / 100.0
    return result


def action(args):
    with named_tempfile('rw') as tfile, args.out as outfile:
        # If query or library file is empty, don't bother executing ssearch.
        # Just print empty file
        # with a header and exit
        if args.fieldnames:
            args_fieldnames = args.fieldnames.split(',')
            fieldnames = []
            for b in USEARCH_HEADER:
                ssearch = toSsearch[b]
                if ssearch == 'qcovs':
                    # will calculate later
                    continue
                elif ssearch in args_fieldnames:
                    fieldnames.append(ssearch)
                else:
                    fieldnames.append(b)
        else:
            fieldnames = USEARCH_HEADER

        if os.stat(args.query).st_size == 0 or \
           os.stat(args.library).st_size == 0:
            # write empty header
            writer = csv.DictWriter(args.out,
                                    extrasaction='ignore',
                                    fieldnames=fieldnames)
            if args.header:
                writer.writeheader()
            return

        command = [args.usearch,
                   '-usearch_global', args.query,
                   '-threads', str(args.threads),
                   '-id', str(args.min_identity),
                   '-db', args.library,
                   '-strand', args.strand,
                   '-blast6out', tfile.name,
                   '-maxaccepts', '0',  # run full algo
                   '-maxhits', args.max or '0']  # '0' = all hits

        log.info(' '.join(command))

        usearch_proc = Popen(command, stderr=PIPE, stdout=PIPE)

        error = set(e.strip() for e in usearch_proc.stderr)
        error = ', '.join(error)

        if usearch_proc.wait() != 0:
            raise CalledProcessError(usearch_proc.returncode, error)

        if error:
            log.error(error)

        usearch_proc.communicate()

        tfile.flush()
        tfile.seek(0)

        results = csv.DictReader(tfile,
                                 fieldnames=fieldnames,
                                 delimiter='\t')

        if 'qcovs' in args_fieldnames or args.min_coverage:
            results = (coverage(r) for r in results)

            if args.min_coverage:
                results = (r for r in results
                           if r['qcovs'] >= args.min_coverage)

        if 'sw_ident' in args_fieldnames:
            # convert to ssearch format
            results = (sw_ident_to_decimal(r) for r in results)

        writer = csv.DictWriter(outfile,
                                fieldnames=args_fieldnames,
                                extrasaction='ignore')

        if args.header:
            writer.writeheader()

        writer.writerows(results)
