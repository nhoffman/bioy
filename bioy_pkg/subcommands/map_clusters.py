"""Create a readmap and clustermap and/or weights file from a
'usearch -uc' cluster file

Input is a usearch6 .uc file

"""

import logging
import sys
import csv
import os

from operator import itemgetter

from bioy_pkg.sequtils import UCLUST_HEADERS, fastalite
from bioy_pkg.utils import Opener, groupbyl

log = logging.getLogger(__name__)


def build_parser(parser):
    parser.add_argument(
        'fastafile', type=lambda f: fastalite(
            Opener()(f), readfile=False),
        help='input fasta file containing original clustered reads')
    parser.add_argument(
        'clusters', type=Opener(),
        help='Clusters file (output of "usearch -uc")')
    parser.add_argument(
        '-o', '--out', type=Opener('w'), default=sys.stdout,
        help='Output fasta containing centroids')
    parser.add_argument(
        '--readmap', type=Opener('w'),
        help='Output file with columns (readname,centroidname)')
    parser.add_argument(
        '--clustermap', type=Opener('w'),
        help='Output file with columns (clustername,samplename)')
    parser.add_argument(
        '--specimen', metavar='SAMPLENAME',
        help='provides samplename for mapfile')
    parser.add_argument(
        '-w', '--weights', type=Opener('w'),
        help='Output file with columns (clustername,weight)')
    parser.add_argument(
        '--min-clust-size', type=int,
        default=1, help='default %(default)s')


def get_mapping(rows):
    """
    target_label (the centroid name) is '*' in rows with type of 'C'
    and 'S'; set the centroid name to query_label (the centroid
    itself) for 'S' rows and drop 'C' rows. This retains clusters of
    size 1 and has the effect of adding the centroid name to the
    list of reads. Returns iterator of (centroid_name, read_name).
    """

    for row in rows:
        if row['type'] == 'C':
            continue
        elif row['type'] == 'S':
            yield (row['query_label'], row['query_label'])
        else:
            yield (row['target_label'], row['query_label'])


def action(args):

    rows = csv.DictReader(args.clusters, delimiter='\t', fieldnames=UCLUST_HEADERS)
    grouped = groupbyl(get_mapping(rows), key=itemgetter(0))  # group by centroid
    clusters = {c: rows for c, rows in grouped if len(rows) >= args.min_clust_size}

    # filter non centroid seqs
    centroids = (c for c in args.fastafile if c.id in clusters)
    for c in centroids:
        args.out.write('>{}\n{}\n'.format(c.description, c.seq))

    readmap = csv.writer(args.readmap) if args.readmap else None
    clustermap = csv.writer(args.clustermap) \
                 if (args.clustermap and args.specimen) else None
    weights = csv.writer(args.weights) if args.weights else None

    for centroid, cluster in clusters.iteritems():
        log.info('writing {}'.format(centroid))

        if readmap:
            readmap.writerows((read, centroid) for read in cluster)

        if clustermap:
            clustermap.writerow((centroid, args.specimen))

        if weights:
            weights.writerow((centroid, str(len(cluster))))

