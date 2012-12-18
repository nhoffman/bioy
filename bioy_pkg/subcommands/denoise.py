"""
Denoise a fasta file of clustered sequences

Input is fasta file containing original reads and the .uc file

Group reads by cluster, run muscle on each, and calculate consensus
sequences.
"""

import logging
import sys
import argparse
import csv
from itertools import groupby, islice
from random import shuffle
from collections import defaultdict

from ion_tools.sequtils import consensus, run_muscle, parse_uc, from_ascii, fastalite
from ion_tools.utils import chunker, opener, Opener

log = logging.getLogger(__name__)

def grouper(groups, rledict = None, min_clust_size = 0, max_clust_size = sys.maxint):
    """
    Return iterator of (seqlist, rlelist) tuples. Clusters are broken
    into chunks no larger than 0.5*max_clust_size
    """

    for cluster_id, cluster in groups:
        # must be a list to determine the cluster size
        cluster = list(cluster)
        if len(cluster) < min_clust_size:
            continue
        elif len(cluster) > max_clust_size:
            shuffle(cluster)
            # combine trailing chunk with next to last if less than
            # half the target chunk size
            combine_last = max_clust_size * 0.5
            for chunk in chunker(cluster, max_clust_size, combine_last):
                rlelist = [rledict[s.id] for s in chunk] if rledict else None
                yield (chunk, rlelist)
        else:
            rlelist = [rledict[s.id] for s in cluster] if rledict else None
            yield (cluster, rlelist)

def build_parser(parser):
    parser.add_argument('fastafile',
            nargs = '?',
            default = sys.stdin,
            help = 'input fasta file containing original clustered reads (default stdin).')
    parser.add_argument('clusters',
            help='Clusters file (output of "usearch -uc")',
            type = Opener())
    parser.add_argument('-r','--rlefile',
            help='An optional file containing run length encoding for infile (.csv.bz2)')
    parser.add_argument('-o','--outfile',
            type=argparse.FileType('w'),
            default = sys.stdout,
            help='Output fasta file.')
    parser.add_argument('-m','--mapfile',
            help='Output file with columns (readname,clustername)')
    parser.add_argument('--max-clust-size',
            type = int,
            default = 100, help = 'default %(default)s')
    parser.add_argument('--min-clust-size',
            type = int,
            default = 10, help = 'default %(default)s')
    parser.add_argument('--limit', metavar = 'N',
                        type = int, help = 'use no more than N seqs')

def action(args):
    log.info('reading %s' % args.clusters)
    cluster_ids, cluster_sizes = parse_uc(args.clusters)

    # sort and group by cluster_id
    get_cluster = lambda s: cluster_ids[s.description]
    seqs = islice(fastalite(opener(args.fastafile)), args.limit)
    groups = groupby(sorted(seqs, key = get_cluster), key = get_cluster)

    log.info('reading %s' % args.rlefile)
    if args.rlefile:
        with opener(args.rlefile) as f:
            rledict = {r['name']:from_ascii(r['rle']) for r in csv.DictReader(f)}
    else:
        rledict = None

    clusters = grouper(groups, rledict,
                       args.min_clust_size, args.max_clust_size)

    # calculate consensus for each cluster, then accumulate names of
    # each set of identical consensus sequences in `exemplars`
    exemplars = defaultdict(list)
    for i, (cluster, rlelist) in enumerate(clusters):
        log.info('aligning cluster %s len %s' % (i, len(cluster)))
        cons = consensus(run_muscle(cluster), rlelist)
        exemplars[cons].extend([s.id for s in cluster])

    # write each consensus sequence
    if args.mapfile:
        mapfile = csv.writer(opener(args.mapfile, 'w'))

    items = sorted(exemplars.items(), key = lambda x: -1*len(x[1]))
    for i, (cons, names) in enumerate(items, start = 1):
        consname = 'cons%04i|%s' % (i, len(names))
        log.info('writing %s' % consname)
        args.outfile.write('>%s\n%s\n' % (consname, cons))

        if args.mapfile:
            mapfile.writerows((name, consname) for name in names)
