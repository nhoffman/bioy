"""
Partition reads in a fastq file by barcode and write an annotated fasta file
"""

import logging
import sys
import csv
import re
from collections import Counter
from difflib import SequenceMatcher
import argparse
import bz2

from Bio import SeqIO

from bioy_pkg.sequtils import homoencode, to_ascii

log = logging.getLogger(__name__)

def get_barcodes(fname):
    """
    Return a list of barcode pairs given a csv file. Assumes
    a header row and that the sample and barcode are in the first and
    second columns respectively.
    """

    with open(fname) as f:
        reader = csv.reader(f)
        reader.next()
        return [(r[1], r[0]) for r in reader]

def find_barcode(seq, bc_rexp, bc_seqs, exact = False, bc_len = 10, search_window = 12):
    """
    Given a SeqRecord, regular expression matching all barcodes, and
    corresponding barcode sequences, return ("barcode","ratio","bc_start","bc_stop")
    """

    seqstr = str(seq.seq)
    substr = seqstr[:search_window]

    min_ratio = (bc_len - 1.0)/bc_len
    min_diff = 3.0/(search_window + bc_len)

    # first look for an exact match to the barcode in the beginning of the
    # sequence.
    match = bc_rexp.search(substr)
    if match:
       bc_start, bc_stop = match.span()
       barcode = match.group()
       ratio = 1.0
    elif exact:
        barcode = ratio = bc_start = bc_stop = None
    else:
        barcode, ratio, bc_start, bc_stop = matcher(substr, bc_seqs, min_ratio, min_diff)

    return barcode, ratio, bc_start, bc_stop

def matcher(seqstr, barcodes, min_ratio = 0, min_diff = 0):
    matches = [SequenceMatcher(a=seqstr, b=barcode) for barcode in barcodes]

    m_ratios = [(m.ratio(), m) for m in matches]
    r0, m0 = max(m_ratios)
    m_ratios.pop(m_ratios.index((r0, m0)))
    r1, m1 = max(m_ratios)

    if r0 - r1 < min_diff:
        return None, None, None, None

    # extract matching substring and calculate ratio
    # Each triple is of the form (i, j, n), and means that a[i:i+n] == b[j:j+n]
    blocks = m0.get_matching_blocks()
    bc_start, bc_stop = blocks[0][1], blocks[-1][1]
    m2 = SequenceMatcher(a=seqstr[bc_start:bc_stop], b=m0.b)
    ratio = m2.ratio()

    if ratio < min_ratio:
        return None, None, None, None

    return m0.b, ratio, bc_start, bc_stop

def build_parser(parser):
    parser.add_argument('fastq',
            help = 'input fastq')
    parser.add_argument('barcodes',
            help = 'csv file containing sample tags (first column) and barcodes (second column); a header row is assumed.')
    parser.add_argument('outfile', help = 'output fasta file',
            type = argparse.FileType('w'))
    parser.add_argument('mapfile',
            help = 'csv.bz2 file with columns "name","label","barcode","ratio","bc_start","bc_stop","rle"')
    parser.add_argument('--stats', help = 'file containing read match statistics',
            type = argparse.FileType('w'), default = sys.stdout)
    parser.add_argument('--unmatched',
            help = 'fasta file containing reads meeting length filters with unmatched barcodes',
            type = argparse.FileType('w'))
    parser.add_argument('--min-length',
            type = int, default = 0,
            help = 'minimum read length')
    parser.add_argument('--max-length',
            type = int,
            default = sys.maxint,
            help = 'maximum read length')
    parser.add_argument('--limit',
            type = int,
            default = sys.maxint,
            help = 'maximum number of query sequences to read from the alignment')
    parser.add_argument('-w', '--search-window',
            type = int,
            default = 12,
            help = 'number of nucleotides from read start to include in barcode scan [%(default)s]')

def action(args):
    barcodes = get_barcodes(args.barcodes)
    bc_seqs, bc_labels = zip(*barcodes)
    bc_dict = dict(barcodes)
    bc_rexp = re.compile('|'.join(bc_seqs))

    opener = bz2.BZ2File if args.fastq.endswith('.bz2') else open
    count = Counter()
    with opener(args.fastq) as f, bz2.BZ2File(args.mapfile, 'w') as mapout:
        writer = csv.writer(mapout)
        writer.writerow(['name','label','barcode','ratio','bc_start','bc_stop','rle'])

        for seq in SeqIO.parse(f, 'fastq'):
            if not args.min_length <= len(seq) <= args.max_length:
                count['fail_len'] += 1
                continue

            barcode, ratio, bc_start, bc_stop = find_barcode(
                seq, bc_rexp, bc_seqs,
                exact = False,
                bc_len = 10,
                search_window = args.search_window)
            sample = bc_dict.get(barcode)

            name = seq.id.replace(':','_')
            if barcode:
                if count['matched'] >= args.limit:
                    break
                count['matched'] += 1
                count[barcode] += 1
                seq, counts = homoencode(seq.seq)
                args.outfile.write('>%s %s\n%s\n' % (name, sample, seq))
                writer.writerow([name, sample, barcode, ratio, bc_start, bc_stop, to_ascii(counts)])
            elif args.unmatched:
                count['no_match'] += 1
                args.unmatched.write('>%s\n%s\n' % (name, seq.seq))
                writer.writerow([name, sample, None, barcode, ratio, bc_start, bc_stop])

    stats = csv.writer(args.stats)
    rows = sorted((bc_dict[k], v) for k, v in count.items() if k in bc_dict)
    stats.writerows(rows)
    stats.writerow(['total_matched', count['matched']])
    stats.writerow(['no_match', count['no_match']])
    stats.writerow(['fail_len', count['fail_len']])
