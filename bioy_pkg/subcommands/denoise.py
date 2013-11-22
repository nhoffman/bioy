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

from operator import itemgetter
from itertools import groupby, islice, chain
from random import shuffle
from collections import defaultdict, Counter
from multiprocessing import Pool

from bioy_pkg.sequtils import consensus, run_muscle, parse_uc, fastalite, from_ascii, homodecode
from bioy_pkg.utils import chunker, Opener, Csv2Dict

log = logging.getLogger(__name__)

CLUSTER_NAME_DELIMITER = '_'

def build_parser(parser):
    parser.add_argument('fastafile',
                        metavar = 'FILE',
                        default = sys.stdin,
                        nargs = '?',
                        type = Opener(),
                        help = 'input fasta file containing original clustered reads (default stdin).')
    parser.add_argument('--clusters',
                        metavar = 'FILE',
                        type = Opener(),
                        help = 'Clusters file (output of "usearch -uc")')
    parser.add_argument('-r','--rlefile', metavar='FILE',
                        type = Csv2Dict('name', 'rle', fieldnames=['name','rle']),
                        help='An optional file containing run length encoding for infile (.csv.bz2)')
    parser.add_argument('-g','--groups', metavar='FILE', type = Opener(),
                        help="""An optional file defining groups for
                             partitioning input reads. If provided,
                             cluster weights will be normalized to
                             proportionally represent each group. File
                             is a headerless csv with columns
                             "seqname","group" (.csv.bz2)""")
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
    parser.add_argument('--max-clust-size', metavar='INTEGER', type = int,
                        default = sys.maxint, help = 'default %(default)s')
    parser.add_argument('--min-clust-size', metavar='INTEGER', type = int,
                        default = 1, help = 'default %(default)s')
    parser.add_argument('--limit', metavar = 'INTEGER', type = int,
                        help = 'use no more than N seqs')
    parser.add_argument('--name-prefix', metavar='STRING',
                        help = 'A string to prepend to each cluster name.')
    parser.add_argument('--name-suffix', metavar='STRING',
                        help = 'A string to append to each cluster name.')
    parser.add_argument('--name-delimiter', default = CLUSTER_NAME_DELIMITER,
                        metavar = 'CHAR',
                        help = 'A character used to delimit elements in cluster names [default "%(default)s"]')

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
        log.debug('aligning cluster {} len {}'.format(i, len(cluster)))
        cons = consensus(run_muscle(cluster), rlelist)

    return cluster, cons


def action(args):

    if args.clusters:
        _, fileExt, = os.path.basename(args.clusters.name).split('.')

        if fileExt == 'uc':
            clusters = parse_uc(args.clusters)[0]
        else:
            clusters = {seq: tag for seq,tag in csv.reader(args.clusters)}

        by_clusters = lambda s: clusters.get(s.id, s.id)
    else:
        by_clusters = lambda _: 'all one cluster'

    seqs = fastalite(args.fastafile)
    seqs = islice(seqs, args.limit)
    seqs = sorted(seqs, key = by_clusters)
    grouped_seqs = groupby(seqs, key = by_clusters)

    chunks = ichunker((group for _, group in grouped_seqs),
                      args.rlefile, args.min_clust_size, args.max_clust_size)

    # calculate consensus for each cluster, then accumulate names of
    # each set of identical consensus sequences in `exemplars` (keys
    # are the consensus sequences themselves).
    exemplars = defaultdict(list)
    pool = Pool(processes = args.threads)
    for cluster, cons in pool.imap_unordered(align_and_consensus, enumerate(chunks, start = 1)):
        exemplars[cons].extend([c.id for c in cluster])

    # calculate ratios of reads for the smallest group to each of the
    # other groups. outseqs is a list of (weight, consensus, list_of_names)
    if args.groups and exemplars:
        groups = dict(csv.reader(args.groups))
        group_counts = Counter(groups[name] for name in chain.from_iterable(exemplars.values()))
        most_common = group_counts.most_common()
        _, least_common = most_common[-1]
        weights = {k: float(least_common)/v for k, v in most_common}
        outseqs = [(sum(weights[groups[n]] for n in names), cons, names)
                   for cons, names in exemplars.items()]
    else:
        outseqs = [(len(names), cons, names) for cons, names in exemplars.items()]

    # write each consensus sequence in descending order of weight
    outseqs.sort(reverse=True, key=itemgetter(0))
    for i, (weight, cons, names) in enumerate(outseqs, start=1):

        name_elements = [args.name_prefix,
                         'cons{:04}'.format(i),
                         '{:.0f}'.format(weight),
                         args.name_suffix]

        consname = args.name_delimiter.join([e for e in name_elements if e])

        log.debug('writing {}'.format(consname))

        args.outfile.write('>{}\n{}\n'.format(consname, cons))

        if args.readmap:
            args.readmap.writerows((name, consname) for name in names)

        if args.clustermap and args.specimen:
            args.clustermap.writerow((consname, args.specimen))

        if args.weights:
            args.weights.writerow((consname, weight))
