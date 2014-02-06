#!/usr/bin/env python

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
    parser.add_argument('infile',
            type = Opener(),
            nargs = '?',
            default = sys.stdin,
            help = """csv file with ssearch36 columns
                      [q_name,q_seq,t_name,t_seq,q_al_start,
                      q_al_stop,t_al_start,t_al_stop,sw_zscore]""")
    parser.add_argument('-i', '--info',
            type = Opener(),
            metavar = 'CSV',
            help = 'info file mapping seqname to tax_id')
    parser.add_argument('-t', '--taxonomy',
            metavar = 'CSV',
            type = Opener(),
            help = 'taxonomy file mapping tax_id to taxonomy')
    parser.add_argument('-o', '--out',
            metavar = 'CSV',
            type = Opener('w'),
            default = sys.stdout,
            help = 'csv output of bases {tax_id, species, positions,,}')
    parser.add_argument('-r', '--rank',
            default = 'species',
            help = 'Aggregate primer stats by specified rank. [%(default)s]')
    parser.add_argument('-f', '--position-freq',
            metavar = 'FLOAT',
            default = 0.05,
            type = float,
            help = 'Minimum base frequency reported for a position [%(default)s]')
    parser.add_argument('-z', '--min-zscore',
            type = float,
            help = 'Minimum z-score value to include alignment in base count.',
            default = 0)

def action(args):
    ### organize seq and tax info
    tax_info = None
    if args.info and args.taxonomy:
        tax = {t['tax_id']:t for t in csv.DictReader(args.taxonomy)}
        tax_info = {i['seqname']:tax[i['tax_id']] for i in csv.DictReader(args.info)}

    ### helper functions
    def intify(al):
        for k in ['q_al_start','q_al_stop','t_sq_len','t_al_start','t_al_stop']:
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

    ### setup and filtering
    aligns = [intify(a) for a in DictReader(args.infile)]

    # assert t_seq uniformity
    aligns = [remove_gaps(a) for a in aligns]
    t_seq = set(a['t_seq'] for a in aligns)
    assert len(t_seq) == 1
    t_seq = t_seq.pop()

    idd_rank = lambda a: get_rank_id(a)
    aligns = sorted(aligns, key = idd_rank)

    total = len(aligns)

    group_by = groupby(aligns, idd_rank)
    total_by_idd = {idd: sum(1 for a in al) for (idd,_), al in group_by}

    # zscore filter
    aligns = [a for a in aligns if a['sw_zscore'] >= args.min_zscore]

    log.info('dropping {} alignments under zscore threshold'.format(total - len(aligns)))

    # alignment within q_seq range
    aligns = [a for a in aligns if in_range(a)]

    log.info('dropping {} partial alignments'.format(total - len(aligns)))
    ###

    ### position count
    bases = defaultdict(Counter)

    for al in aligns:
        idd, rank = get_rank_id(al)

        # reference sequence, starting at first base of primer alignment
        q_start = al['q_al_start'] - 1
        q_stop = al['q_al_stop']
        qseq = al['q_seq'][q_start:q_stop]
        for pos, base in enumerate(qseq):
            bases[(pos, idd, rank)][base] += 1

    fieldnames = ['name'] if tax_info else []
    fieldnames += ['position', 'A', 'T', 'G', 'C', 'N',
                   'expected', 'alignments', 'total']
    fieldnames += ['rank', 'id'] if tax_info else []

    out = DictWriter(args.out, fieldnames = fieldnames, extrasaction='ignore')
    out.writeheader()

    organized_bases = OrderedDict(sorted(bases.items(), key = lambda b: (b[0][1], b[0][0])))

    for (pos, idd, rank), counter in organized_bases.items():
        naligns = sum([v for v in counter.values()])
        out.writerow({
            'name':tax[idd]['tax_name'] if tax_info else '',
            'rank':rank,
            'id':idd,
            'position':pos+1,
            'expected':t_seq[pos],
            'total':total_by_idd[idd],
            'alignments':naligns,
            'A':pop_base_or_zero(counter, 'A'),
            'T':pop_base_or_zero(counter, 'T'),
            'G':pop_base_or_zero(counter, 'G'),
            'C':pop_base_or_zero(counter, 'C'),
            'N':sum(counter.values())
            })

