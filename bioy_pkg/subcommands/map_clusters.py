"""
Create a readmap and clustermap and/or weights file from a 'usearch -uc' cluster file

Input is a usearch6 .uc file
"""

import logging
import sys
import csv

from itertools import groupby
from operator import itemgetter

from bioy_pkg.sequtils import UCLUST_HEADERS, fastalite
from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('--fasta',
            type = lambda f: fastalite(Opener()(f), readfile = False),
            help = 'input fasta file containing original clustered reads')
    parser.add_argument('clusters',
            type = Opener(),
            help = 'Clusters file (output of "usearch -uc")')
    parser.add_argument('--specimen',
            help = 'sample name for specimen map')
    parser.add_argument('--fasta-out',
            type = Opener('w'),
            help = 'Output fasta mapping reads to centroid (readname,centroidname)')
    parser.add_argument('--out',
            type = lambda f: csv.writer(Opener('w')(f)),
            default = sys.stdout,
            help = 'Output file with columns (readname,samplename)')
    parser.add_argument('--specimen-map',
            type = lambda f: csv.writer(Opener('w')(f)),
            help = 'Output file with columns (clustername,samplename)')
    parser.add_argument('-w', '--weights',
            type = lambda f: csv.writer(Opener('w')(f)),
            help = 'Output file with columns (clustername,weight)')
    parser.add_argument('--min-clust-size',
            type = int,
            default = 1, help = 'default %(default)s')

def action(args):
    clusters = csv.DictReader(args.clusters, delimiter = '\t', fieldnames = UCLUST_HEADERS)
    clusters = sorted(clusters, key = itemgetter('target_label'))
    clusters = groupby(clusters, key = itemgetter('target_label'))
    clusters = ((c,list(r)) for c,r in clusters)
    clusters = ((c,r) for c,r in clusters if len(r) >= args.min_clust_size)
    clusters = dict(clusters)
    clusters.pop('*') # remove star seqs (centroids are the rest of keys)

    for centroid, cluster in clusters.items():
        log.info('writing {}'.format(centroid))

        args.out.writerows((data['query_label'], centroid) for data in cluster)

        if args.specimen_map and args.specimen:
            args.specimen_map.writerow((centroid, args.specimen))

        if args.weights:
            weight = str(len(cluster))
            args.weights.writerow((centroid, weight))

    # filter non centroid seqs
    if args.fasta_out:
        for c in  (c for c in args.fasta if c.id in clusters):
            args.fasta_out.write('>{}\n{}\n'.format(c.description, c.seq))


