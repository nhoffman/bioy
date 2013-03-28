"""
Prepare a table annotating specimens for the run
"""

import logging
import sys
import csv

from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument('labels',
                        default = sys.stdin, type = Opener('rU'),
                        help = ('csv format file (minimally) with '
                                'headers "barcode_id,label"'))
    parser.add_argument('-o', '--out',
                        type = Opener('w'),
                        default = sys.stdout,
                        help = 'csv file with headers "specimen,label,..."')
    parser.add_argument('run',
                        help = 'string identifying run in specimen labels')


def action(args):

    reader = csv.DictReader(args.labels)
    fields_in = [f.lower() for f in reader.fieldnames]
    reader.fieldnames = fields_in

    required = ['barcode_id', 'label']
    if not set(required).issubset(set(fields_in)):
        sys.exit('the following filednames are required: '.format(
            ''.join(required)))

    fields_out = ['specimen'] + fields_in

    writer = csv.DictWriter(args.out, fieldnames = fields_out)
    writer.writeheader()
    for d in reader:
        d['specimen'] = '{}s{:02}'.format(args.run, int(d['barcode_id']))
        writer.writerow(d)
