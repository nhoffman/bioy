"""
Return the children of a taxtable given a list of taxids
"""

import csv
import logging
import sys

from operator import itemgetter

from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('taxids',
            type = lambda l: set(l.split(',')),
            help = 'comma delimited list of column names to group columns')
    parser.add_argument('taxonomy',
            metavar = 'CSV',
            default = sys.stdin,
            nargs = '?',
            type = Opener(),
            help = 'input classify csv file')
    parser.add_argument('-o', '--out',
            type = Opener('w'),
            default = sys.stdout,
            help = """csv with columns
            [parent_name,parent_id,parent_rank,tax_name,tax_id,rank]""")

def action(args):
    taxonomy = csv.DictReader(args.taxonomy)
    ranks = taxonomy.fieldnames[4:]

    tax = {}

    # filter out claves
    dropped = total = 0
    for t in taxonomy:
        parent_id = set(itemgetter(*ranks)(t)) & args.taxids
        if parent_id:
            t['parent_id'], = parent_id
            tax[t['tax_id']] = t
        else:
            dropped += 1
        total += 1

    log.info('dropped {} of {} records ({:.2%}) as not children'.format(
        dropped, total, float(dropped) / total))

    # add parent_name and parent_rank
    for t in tax.values():
        parent = tax[t['parent_id']]
        t['parent_name'] = parent['tax_name']
        t['parent_rank'] = parent['rank']

    # sort by parent_name and tax_name ??
    rows = sorted(tax.values(), key = itemgetter('parent_name', 'tax_name'))

    # output
    fieldnames = ['parent_name', 'parent_id', 'parent_rank',
                  'tax_name', 'tax_id', 'rank']
    out = csv.DictWriter(args.out,
                         extrasaction = 'ignore',
                         fieldnames = fieldnames)
    out.writeheader()
    out.writerows(rows)

