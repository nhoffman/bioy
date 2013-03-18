"""
Outputs a standard Genbank Record File into fasta file format and optional seqinfo file in format ['seqname','tax_id','accession','description','length','ambig_count','is_type','rdp_lineage']
"""

import sys

from bioy_pkg.utils import Opener
from Bio import SeqIO
from csv import DictWriter

from itertools import ifilter, islice

from bioy_pkg.sequtils import UNCLASSIFIED_REGEX, INFO_HEADER, is_type, gb2info

def build_parser(parser):
    parser.add_argument('infile',
            nargs = '?',
            default = sys.stdin,
            type = Opener('r'),
            help = 'fasta file')
    parser.add_argument('-o', '--out',
            default = sys.stdout,
            type = Opener('w'),
            help = 'ouput fasta file')
    parser.add_argument('-O', '--info-out',
            type = lambda f: DictWriter(Opener('w')(f), fieldnames = INFO_HEADER),
            help = 'Output seq info file')
    parser.add_argument('-f', '--filter',
            action = 'store_true',
            help = 'Filter against UNCLASSIFIED_REGEX: {}'.format(UNCLASSIFIED_REGEX))
    parser.add_argument('--products',
            action = 'append',
            help = 'parse ribosomal subunits from genbank record')
    parser.add_argument('-t', '--type-strains',
            action = 'store_true',
            help = 'Only return type strain sequences')
    parser.add_argument('-l', '--limit',
            type = int,
            help = 'limit the number of records to load')
    parser.add_argument('-m', '--minus',
            action = 'store_true',
            help = 'seq is in minus strand orientation')

def action(args):
    records = islice(SeqIO.parse(args.infile, 'genbank'), args.limit)

    if args.type_strains:
        records = ifilter(lambda r: is_type(r), records)

    if args.filter:
        records = ifilter(lambda r: not UNCLASSIFIED_REGEX.search(r.description), records)

    info = []

    if args.products:
        for r in records:
            for f in r.features:
                products = f.qualifiers.get('product', [])
                if any(r in products for r in args.products):
                    tag = f.qualifiers.get('locus_tag', ['unspecified'])[0]
                    start, end = f.location.start.position, f.location.end.position
                    seq = r.seq[start:end]
                    name = '{}_{}_{}'.format(tag, start, end)

                    if (args.minus and f.location.strand == 1) or f.location.strand == -1:
                        seq = seq.reverse_complement()

                    args.out.write('>{} {} {}\n{}\n'.format(name, r.id, r.description, seq))
                    info.append(gb2info(r, seqname = name))
    else:
        for r in records:
            src = next(ifilter(lambda f: f.type == 'source', r.features), False)
            if src:
                seq = r.seq

                if (args.minus and src.location.strand == 1) or src.location.strand == -1:
                    seq = seq.reverse_complement()

                args.out.write('>{} {} {}\n{}\n'.format(r.name, r.id, r.description, seq))
                info.append(gb2info(r, seqname = r.name))

    if args.info_out:
        args.info_out.writeheader()
        args.info_out.writerows(info)

