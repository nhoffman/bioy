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
Tally ssearch base count by position

Groups results by taxonomy rank
"""

import csv
import logging
import sys

from collections import defaultdict, OrderedDict, Counter
from csv import DictReader, DictWriter
from itertools import groupby

from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument(
        'infile',
        type=Opener(),
        nargs='?',
        default=sys.stdin,
        help=('csv file with ssearch36 columns '
              '[q_name,q_seq,t_name,t_seq,q_al_start,q_al_stop,'
              't_al_start,t_al_stop,t_sq_len,sw_zscore]'))
    parser.add_argument(
        '-i', '--info',
        type=Opener(),
        metavar='CSV',
        help='info file mapping seqname to tax_id')
    parser.add_argument(
        '-t', '--taxonomy',
        metavar='CSV',
        type=Opener(),
        help='taxonomy file mapping tax_id to taxonomy')
    parser.add_argument(
        '-o', '--out',
        metavar='CSV',
        type=Opener('w'),
        default=sys.stdout,
        help='csv output of bases {tax_id, species, positions,,}')
    parser.add_argument(
        '-r', '--rank',
        default='species',
        help='Aggregate primer stats by specified rank. [%(default)s]')
    parser.add_argument(
        '-f', '--position-freq',
        metavar='FLOAT',
        default=0.05,
        type=float,
        help='Minimum base frequency reported for a position [%(default)s]')
    parser.add_argument(
        '-z', '--min-zscore',
        type=float,
        help='Minimum z-score value to include alignment in base count.',
        default=0)


def action(args):
    # organize seq and tax info
    tax_info = None
    if args.info and args.taxonomy:
        tax = {t['tax_id']: t for t in csv.DictReader(args.taxonomy)}
        tax_info = {i['seqname']: tax[i['tax_id']]
                    for i in csv.DictReader(args.info)}

    # helper functions
    def intify(al):
        for k in ['q_al_start', 'q_al_stop',
                  't_sq_len', 't_al_start', 't_al_stop']:
            al[k] = int(al[k])
        return al

    def get_rank_id(a):
        rank = idd = None

        if tax_info:
            idd = tax_info[a['q_name']][args.rank]
            if idd:
                rank = args.rank
            else:
                rank = 'species'
                idd = tax_info[a['q_name']][rank]

        return idd, rank

    def in_range(a):
        return a['q_al_stop'] - a['q_al_start'] == a['t_sq_len'] - 1

    def pop_base_or_zero(counter, base):
        return counter.pop(base) if base in counter else 0

    def remove_gaps(a):
        a['t_seq'] = a['t_seq'].replace('-', '')
        return a

    # setup and filtering
    aligns = [intify(a) for a in DictReader(args.infile)]

    # assert t_seq uniformity
    aligns = [remove_gaps(a) for a in aligns]
    t_seq = set(a['t_seq'] for a in aligns)
    assert len(t_seq) == 1
    t_seq = t_seq.pop()

    aligns = sorted(aligns, key=get_rank_id)

    total = len(aligns)

    group_by = groupby(aligns, get_rank_id)
    total_by_idd = {idd: sum(1 for a in al) for (idd, _), al in group_by}

    # zscore filter
    aligns = [a for a in aligns if float(a['sw_zscore']) >= args.min_zscore]

    log.info(
        'dropping {} alignments under zscore threshold'.format(
            total - len(aligns)))

    # alignment within q_seq range
    aligns = [a for a in aligns if in_range(a)]

    log.info('dropping {} partial alignments'.format(total - len(aligns)))
    ###

    # position count
    bases = defaultdict(Counter)

    for al in aligns:
        idd, rank = get_rank_id(al)

        # reference sequence, starting at first base of primer alignment
        q_start = al['q_al_start'] - 1
        q_stop = al['q_al_stop']
        qseq = al['q_seq'][q_start:q_stop]
        for pos, base in enumerate(qseq):
            bases[(pos, idd, rank)][base] += 1

    fieldnames = ['tax_name'] if tax_info else []
    fieldnames += ['position', 'A', 'T', 'G', 'C', 'N',
                   'expected', 'naligns', 'nseqs']
    fieldnames += ['rank', 'tax_id'] if tax_info else []

    out = DictWriter(args.out, fieldnames=fieldnames, extrasaction='ignore')
    out.writeheader()

    organized_bases = OrderedDict(
        sorted(
            bases.items(),
            key=lambda b: (
                b[0][1],
                b[0][0])))

    for (pos, idd, rank), counter in organized_bases.items():
        naligns = sum([v for v in counter.values()])
        out.writerow({
            'tax_name': tax[idd]['tax_name'] if tax_info else '',
            'rank': rank,
            'tax_id': idd,
            'position': pos + 1,
            'expected': t_seq[pos],
            'nseqs': total_by_idd[idd],
            'naligns': naligns,
            'A': pop_base_or_zero(counter, 'A'),
            'T': pop_base_or_zero(counter, 'T'),
            'G': pop_base_or_zero(counter, 'G'),
            'C': pop_base_or_zero(counter, 'C'),
            'N': sum(counter.values())
        })
