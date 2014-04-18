"""Create a readmap and clustermap and/or weights file from a
'usearch -uc' cluster file

Input is a usearch6 .uc file

"""

import logging
import sys
import csv
import os

from itertools import groupby
from operator import itemgetter

from bioy_pkg.sequtils import UCLUST_HEADERS, fastalite
from bioy_pkg.utils import Opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('fastafile',
            type = lambda f: fastalite(Opener()(f), readfile = False),
            help = 'input fasta file containing original clustered reads')
    parser.add_argument('clusters', type = Opener(),
            help = 'Clusters file (output of "usearch -uc")')
    parser.add_argument('-o','--out', type = Opener('w'), default = sys.stdout,
            help = 'Output fasta containing centroids')
    parser.add_argument('--readmap', type = Opener('w'),
            help = 'Output file with columns (readname,centroidname)')
    parser.add_argument('--clustermap', type = Opener('w'),
            help = 'Output file with columns (clustername,samplename)')
    parser.add_argument('--specimen', metavar='SAMPLENAME',
            help = 'provides samplename for mapfile')
    parser.add_argument('-w', '--weights', type = Opener('w'),
            help = 'Output file with columns (clustername,weight)')
    parser.add_argument('--min-clust-size', type = int,
            default = 1, help = 'default %(default)s')

def action(args):
    clusters = csv.DictReader(args.clusters, delimiter = '\t', fieldnames = UCLUST_HEADERS)
    clusters = sorted(clusters, key = itemgetter('target_label'))
    clusters = groupby(clusters, key = itemgetter('target_label'))
    clusters = ((c,list(r)) for c,r in clusters)
    clusters = ((c,r) for c,r in clusters if len(r) >= args.min_clust_size)
    clusters = dict(clusters)
    clusters.pop('*') # remove star seqs (centroids are the rest of keys)

    # filter non centroid seqs
    centroids = (c for c in args.fastafile if c.id in clusters)

    for c in centroids:
        args.out.write('>{}\n{}\n'.format(c.description, c.seq))

    readmap = csv.writer(args.readmap) if args.readmap else None
    clustermap = csv.writer(args.clustermap) \
                 if args.clustermap and args.specimen else None
    weights = csv.writer(args.weights) if args.weights else None

    for centroid, cluster in clusters.items():
        log.info('writing {}'.format(centroid))

        if readmap:
            readmap.writerows((data['query_label'], centroid) for data in cluster)

        if clustermap:
            writerow((centroid, args.specimen))

        if weights:
            weights.writerow((centroid, str(len(cluster))))

