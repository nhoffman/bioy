"""
Create a readmap and clustermap and/or weights file from a 'usearch -uc' cluster file

Input is a usearch 6 .uc file
"""

import logging
import sys
import csv
import os

from itertools import groupby
from operator import itemgetter

from bioy_pkg.sequtils import UCLUST_HEADERS
from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('clusters',
            nargs = '?',
            type = Opener(),
            help = 'Clusters file (output of "usearch -uc")')
    parser.add_argument('--specimen',
            help = 'sample name for mapfile')
    parser.add_argument('-o','--out',
            type = Opener('w'),
            default = sys.stdout,
            help = 'Output fasta mapping reads to centroid (readname,centroidname)')
    parser.add_argument('--clustermap',
            type = lambda f: csv.writer(Opener('w')(f)),
            help = 'Output file with columns (clustername,samplename)')
    parser.add_argument('-w', '--weights',
            type = lambda f: csv.writer(Opener('w')(f)),
            help = 'Output file with columns (clustername,weight)')
    parser.add_argument('--min-clust-size',
            type = int,
            default = 1, help = 'default %(default)s')
    parser.add_argument('--threads',
            default = int(os.environ.get('THREADS_ALLOC') or 1),
            type = int,
            help = """Number of threads (CPUs) to use.
                   Can also specify with environment variable THREADS_ALLOC
                   default = %(default)s""")

def action(args):
    clusters = csv.DictReader(args.clusters, delimiter = '\t', fieldnames = UCLUST_HEADERS)
    clusters = sorted(clusters, key = itemgetter('target_label'))
    clusters = groupby(clusters, key = itemgetter('target_label'))
    clusters = ((c,list(r)) for c,r in clusters)
    clusters = ((c,r) for c,r in clusters if len(r) >= args.min_clust_size)
    clusters = ((c,r) for c,r in clusters if c != '*')

    out = csv.writer(args.out)

    for centroid, cluster in clusters:
        log.info('writing {}'.format(centroid))

        out.writerows((data['query_label'], centroid) for data in cluster)

        if args.clustermap and args.specimen:
            args.clustermap.writerow((centroid, args.specimen))

        if args.weights:
            weight = str(len(cluster))
            args.weights.writerow((centroid, weight))
