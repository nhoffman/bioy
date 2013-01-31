"""
Ouput primer alignment stats
"""

import logging
import sys
import csv
from itertools import islice, groupby

from bioy_pkg.sequtils import parse_ssearch36, parse_primer_alignments
from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('alignments',
            type = Opener(),
            help = 'ssearch output (bz2 compressed)')
    parser.add_argument('-o','--outfile',
            help = 'output csv file (appends to existing file)',
            type = Opener('w'),
            default = sys.stdout)
    parser.add_argument('--limit',
            help = 'Quit after N query sequences [%(default)s]',
            type = int)
def action(args):
    writer = csv.DictWriter(
        args.outfile,
        fieldnames = ['name','primer','start','stop','sw_expect',
            'sw_ident','sw_overlap', 'sw_zscore'],
        extrasaction = 'ignore')
    writer.writeheader()

    for query, hits in islice(groupby(parse_ssearch36(args.alignments),
        lambda x: x['q_name']), args.limit):
        aligndata = parse_primer_alignments(hits)
        for primer in ['l','r']:
            d = dict(name=aligndata['name'], primer=primer, **aligndata[primer])
            writer.writerow(d)
