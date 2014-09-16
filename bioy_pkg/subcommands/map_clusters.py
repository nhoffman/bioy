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

"""Create a readmap and specimenmap and/or weights file from a
'usearch -uc' cluster file

Input is a usearch6 .uc file

"""

import logging
import sys
import csv

from collections import Counter
from operator import itemgetter

from bioy_pkg.sequtils import UCLUST_HEADERS, fastalite
from bioy_pkg.utils import Opener, groupbyl

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument(
        'clusters', type=Opener(),
        help='Clusters file (output of "usearch -uc")')
    parser.add_argument(
        '--fasta-in', type=lambda f: fastalite(Opener()(f)),
        help='input fasta file containing original clustered reads')
    parser.add_argument(
        '--fasta-out', type=Opener('w'),
        help='Output fasta containing centroids')
    parser.add_argument(
        '-g', '--groups', metavar='FILE', type=Opener(),
        help="""An optional file defining groups for partitioning
        input reads. If provided, cluster weights will be normalized
        to proportionally represent each group. File is a headerless
        csv with columns "seqname","group" (.csv.bz2)""")
    parser.add_argument(
        '--min-clust-size', type=int,
        default=1, help='[%(default)s]')
    parser.add_argument(
        '-o', '--out', type=Opener('w'), default=sys.stdout,
        help='Output file with columns (readname,centroidname)')
    parser.add_argument(
        '--specimenmap', type=Opener('w'),
        help='Output file with columns (clustername,samplename)')
    parser.add_argument(
        '--specimen', metavar='SAMPLENAME',
        help='provides samplename for mapfile')
    parser.add_argument(
        '-w', '--weights', type=Opener('w'),
        help='Output file with columns (clustername,weight)')

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
    readmap = csv.writer(args.out)
    specimenmap = csv.writer(args.specimenmap) \
                 if (args.specimenmap and args.specimen) else None
    weights = csv.writer(args.weights) if args.weights else None

    # Calculate ratios of reads for the smallest group to each of the
    # other groups.
    if args.groups:
        groups = dict(csv.reader(args.groups))
        most_common = Counter(groups.values()).most_common()
        _, least_common = most_common[-1]
        wdict = {k: float(least_common) / v for k, v in most_common}

    for centroid, cluster in clusters.iteritems():
        log.info('writing {}'.format(centroid))

        for _centroid, read in cluster:
            readmap.writerow([read, centroid])

        if specimenmap:
            specimenmap.writerow((centroid, args.specimen))

        if weights:
            if args.groups:
                # normalize weight of each cluster by contribution of
                # reads by each group defined in --groups
                weights.writerow((centroid, sum(wdict[groups[r]] for c, r in cluster)))
            else:
                weights.writerow((centroid, len(cluster)))

    # filter non centroid seqs
    if args.fasta_in and args.fasta_out:
        for c in  (c for c in args.fasta_in if c.id in clusters):
            args.fasta_out.write('>{}\n{}\n'.format(c.description, c.seq))

