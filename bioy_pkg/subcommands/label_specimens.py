"""
Prepare a table annotating specimens for the run
"""

import logging
import sys
import csv

from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)

def format_label(data_id, barcode_id):
    return '{:03}_{:02}'.format(int(data_id), int(barcode_id))

def build_parser(parser):
    parser.add_argument('labels',
                        default = sys.stdin, type = Opener('rU'),
                        help = ('csv format file (minimally) with '
                                'headers "barcode_id,label"'))
    parser.add_argument('data_id', type = int,
                        help = 'integer identifying the data set')
    parser.add_argument('-o', '--out',
                        type = Opener('w'),
                        default = sys.stdout,
                        help = 'csv file with headers "specimen,label,..."')


def action(args):

    reader = csv.DictReader(args.labels)
    fields_in = [f.lower() for f in reader.fieldnames]
    reader.fieldnames = fields_in

    required = ['barcode_id', 'label']
    if not set(required).issubset(set(fields_in)):
        sys.exit('the following filednames are required: '.format(
            ''.join(required)))

    # remove empty columns
    fields_out = [f for f in ['specimen'] + fields_in if f.strip()]

    writer = csv.DictWriter(args.out, fieldnames = fields_out,
                            extrasaction = 'ignore')
    writer.writeheader()
    for d in reader:
        if d['barcode_id']:
            d['specimen'] = format_label(args.data_id, d['barcode_id'])
            writer.writerow(d)
