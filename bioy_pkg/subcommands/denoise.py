"""
Denoise a fasta file of clustered sequences

Input is fasta file containing original reads and the .uc file

Group reads by cluster, run muscle on each, and calculate consensus
sequences.
"""

import logging
import sys
import csv
import os

from itertools import groupby, islice
from random import shuffle
from collections import defaultdict
from multiprocessing import Pool

from bioy_pkg.sequtils import consensus, run_muscle, parse_uc, fastalite, from_ascii, homodecode
from bioy_pkg.utils import chunker, Opener, Csv2Dict

log = logging.getLogger(__name__)

CLUSTER_NAME_DELIMITER = '_'

def build_parser(parser):
    parser.add_argument('fastafile',
            type = lambda f: fastalite(Opener()(f), readfile = False),
            help = 'input fasta file containing original clustered reads (default stdin).')
    parser.add_argument('clusters',
            type = Opener(),
            help = 'Clusters file (output of "usearch -uc")')
    parser.add_argument('-r','--rlefile', metavar='FILE',
            type = Csv2Dict('name', 'rle', fieldnames=['name','rle']),
            help='An optional file containing run length encoding for infile (.csv.bz2)')
    parser.add_argument('-o','--outfile', metavar='FILE',
            type = Opener('w'),
            default = sys.stdout,
            help='Output fasta file.')
    parser.add_argument('--readmap', metavar='FILE',
            type = lambda f: csv.writer(Opener('w')(f)),
            help = 'Output file with columns (readname,clustername)')
    parser.add_argument('--clustermap', metavar='FILE',
            type = lambda f: csv.writer(Opener('w')(f)),
                        help = 'Output file with columns (clustername,specimen)')
    parser.add_argument('--specimen', metavar='NAME',
            help = 'value for second column of --clustermap')
    parser.add_argument('-w', '--weights', metavar='FILE',
            type = lambda f: csv.writer(Opener('w')(f)),
            help = 'Output file with columns (clustername,weight)')
    parser.add_argument('--max-clust-size', metavar='INTEGER',
            type = int,
            default = sys.maxint, help = 'default %(default)s')
    parser.add_argument('--min-clust-size', metavar='INTEGER',
            type = int,
            default = 1, help = 'default %(default)s')
    parser.add_argument('--limit',
            metavar = 'INTEGER',
            type = int,
            help = 'use no more than N seqs')
    parser.add_argument('--name-prefix', metavar='STRING',
            help = 'A string to prepend to each cluster name.')
    parser.add_argument('--name-suffix', metavar='STRING',
            help = 'A string to append to each cluster name.')
    parser.add_argument('--name-delimiter', default = CLUSTER_NAME_DELIMITER,
                        metavar = 'CHAR',
                        help = 'A character used to delimit elements in cluster names [default "%(default)s"]')
    parser.add_argument('--threads',
            default = int(os.environ.get('THREADS_ALLOC') or 1),
            type = int,
            help = """Number of threads (CPUs) to use.
                   Can also specify with environment variable THREADS_ALLOC
                   default = %(default)s""")

def ichunker(seqs, rledict=None, min_clust_size=1, max_clust_size=sys.maxint):
    """Return iterator of (seqlist, rlelist) tuples. Clusters are broken
    into chunks no larger than 0.5*max_clust_size. Returns an iterator
    of (cluster_seqs, rle_info). If rle_info is not provided, rle_info
    is None.

    """

    for cluster in seqs:
        cluster = list(cluster)
        if len(cluster) < min_clust_size:
            continue
        elif len(cluster) > max_clust_size:
            shuffle(cluster)
            # combine trailing chunk with next to last if less than
            # half the target chunk size
            combine_last = max_clust_size * 0.5
            for chunk in chunker(cluster, max_clust_size, combine_last):
                rlelist = [from_ascii(rledict[s.id]) for s in chunk] if rledict else None
                yield (chunk, rlelist)
        else:
            rlelist = [from_ascii(rledict[s.id]) for s in cluster] if rledict else None
            yield (cluster, rlelist)


def align_and_consensus(chunk):
    """Wraps functions for alignment and consensus generation. Does not
    perform alignment for clusters of length 1.

    """

    i, (cluster, rlelist) = chunk

    if len(cluster) == 1:
        # no need to align...
        seq = cluster[0]
        rle = rlelist[0] if rlelist else None
        cons = homodecode(seq, rle) if rle else seq.seq
    else:
        log.info('aligning cluster {} len {}'.format(i, len(cluster)))
        cons = consensus(run_muscle(cluster), rlelist)

    return cluster, cons


def action(args):

    _, fileExt, = os.path.basename(args.clusters.name).split('.')

    if fileExt == 'uc':
        clusters = parse_uc(args.clusters)[0]
    else:
        clusters = {seq: tag for seq,tag in csv.reader(args.clusters)}

    by_clusters = lambda s: clusters.get(s.id, s.id)

    seqs = islice(args.fastafile, args.limit)
    seqs = sorted(seqs, key = by_clusters)
    groups = groupby(seqs, key = by_clusters)
    chunks = ichunker((group for _, group in groups),
                      args.rlefile, args.min_clust_size, args.max_clust_size)
    chunks = list(chunks)

    # calculate consensus for each cluster, then accumulate names of
    # each set of identical consensus sequences in `exemplars`
    exemplars = defaultdict(list)

    pool = Pool(processes = args.threads)

    for cluster, cons in pool.imap_unordered(align_and_consensus, enumerate(chunks, start = 1)):
        exemplars[cons].extend([c.id for c in cluster])

    # write each consensus sequence in descending order of cluster size
    items = sorted(exemplars.items(), key=lambda x: len(x[1]), reverse=True)
    for i, (cons, names) in enumerate(items, start=1):
        name_elements = [args.name_prefix,
                         'cons{:04}'.format(i),
                         str(len(names)),
                         args.name_suffix]
        consname = args.name_delimiter.join([e for e in name_elements if e])

        log.info('writing {}'.format(consname))

        args.outfile.write('>{}\n{}\n'.format(consname, cons))

        if args.readmap:
            args.readmap.writerows((name, consname) for name in names)

        if args.clustermap and args.specimen:
            args.clustermap.writerow((consname, args.specimen))

        if args.weights:
            args.weights.writerow((consname, weight))
