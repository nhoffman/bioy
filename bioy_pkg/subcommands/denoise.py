"""
Denoise a fasta file of clustered sequences

Input is fasta file containing original reads and the .uc file

Group reads by cluster, run muscle on each, and calculate consensus
sequences.
"""

import logging
import sys
import csv
from itertools import groupby, islice, imap, ifilter
from random import shuffle
from collections import defaultdict

from bioy_pkg.sequtils import consensus, run_muscle, parse_uc, fastalite, from_ascii
from bioy_pkg.utils import chunker, Opener, Csv2Dict

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('fastafile',
            nargs = '?',
            default = sys.stdin,
            type = lambda f: fastalite(Opener()(f)),
            help = 'input fasta file containing original clustered reads (default stdin).')
    parser.add_argument('clusters',
            type = lambda c: parse_uc(Opener()(c))[0],
            help = 'Clusters file (output of "usearch -uc")')
    parser.add_argument('--specimen',
            help = 'sample name for mapfile')
    parser.add_argument('-r','--rlefile',
            type = Csv2Dict(fieldnames=['name','rle']),
            help='An optional file containing run length encoding for infile (.csv.bz2)')
    parser.add_argument('-o','--outfile',
            type = Opener('w'),
            default = sys.stdout,
            help='Output fasta file.')
    parser.add_argument('--readmap',
            type = lambda f: csv.writer(Opener('w')(f)),
            help = 'Output file with columns (readname,clustername)')
    parser.add_argument('--clustermap',
            type = lambda f: csv.writer(Opener('w')(f)),
            help = 'Output file with columns (clustername,samplename)')
    parser.add_argument('-w', '--weights',
            type = lambda f: csv.writer(Opener('w')(f)),
            help = 'Output file with columns (clustername,weight)')
    parser.add_argument('--max-clust-size',
            type = int,
            default = 100, help = 'default %(default)s')
    parser.add_argument('--min-clust-size',
            type = int,
            default = 10, help = 'default %(default)s')
    parser.add_argument('--limit', metavar = 'N',
            type = int,
            help = 'use no more than N seqs')

def ichunker(seqs, rledict = None, max_clust_size = sys.maxint):
    """
    Return iterator of (seqlist, rlelist) tuples. Clusters are broken
    into chunks no larger than 0.5*max_clust_size
    """

    for cluster in seqs:
        if len(cluster) > max_clust_size:
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

def action(args):
    seqs = islice(args.fastafile, args.limit)
    seqs = sorted(seqs, key = lambda s: args.clusters[s.description])
    seqs = groupby(seqs, lambda s: args.clusters[s.description])
    seqs = imap(lambda (_,s): list(s), seqs)
    seqs = ifilter(lambda s: len(s) >= args.min_clust_size, seqs)
    seqs = ichunker(seqs, args.rlefile, args.max_clust_size)

    # calculate consensus for each cluster, then accumulate names of
    # each set of identical consensus sequences in `exemplars`
    exemplars = defaultdict(list)
    for i, (cluster, rlelist) in enumerate(seqs):
        log.info('aligning cluster {} len {}'.format(i, len(cluster)))
        cons = consensus(run_muscle(cluster), rlelist)
        exemplars[cons].extend([s.id for s in cluster])

    # write each consensus sequence
    items = sorted(exemplars.items(), key = lambda x: -1 * len(x[1]))
    for i, (cons, names) in enumerate(items, start = 1):
        weight = len(names)
        consname = 'cons{:04}|{}'.format(i, weight)

        log.info('writing {}'.format(consname))

        args.outfile.write('>{}\n{}\n'.format(consname, cons))

        if args.readmap:
            args.readmap.writerows((name, consname) for name in names)

        if args.clustermap and args.specimen:
            args.clustermap.writerow((consname, args.specimen))

        if args.weights:
            args.weights.writerow((consname, weight))
