"""
Outputs a standard Genbank Record File into fasta file format and optional seqinfo file in format ['seqname','tax_id','accession','description','length','ambig_count','is_type','rdp_lineage']
"""

import sys
import logging

from bioy_pkg.utils import Opener
from Bio import SeqIO
from csv import DictWriter

from itertools import ifilter, islice

from bioy_pkg.sequtils import UNCLASSIFIED_REGEX,tax_of_genbank, count_ambiguous, is_type

log = logging.getLogger(__name__)

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
            type = lambda f: DictWriter(Opener('w')(f), fieldnames =
                ['seqname','tax_id','accession','description','length','ambig_count']),
            help = 'Output seq info file')
    parser.add_argument('-f', '--filter',
            action = 'store_true',
            help = 'Filter against UNCLASSIFIED_REGEX: {}'.format(UNCLASSIFIED_REGEX))
    parser.add_argument('--feature',
            action = 'append',
            dest = 'features',
            help = """parse genome features (ex ribosomal subunits (16S, 18S, etc))
                      from genbank record""")
    parser.add_argument('-t', '--type-strains',
            action = 'store_true',
            help = 'Only return type strain sequences')
    parser.add_argument('-l', '--limit',
            type = int,
            help = 'limit the number of records to load')
    parser.add_argument('-m', '--minus',
            action = 'store_true',
            help = 'seq is in minus strand orientation')
    parser.add_argument('--max-ambiguous',
            type = int,
            default = sys.maxint,
            help = 'ambiguous base count threshold default = %(default)s')
    parser.add_argument('--min-length',
            type = int,
            default = 0,
            help = 'minimum length for reference sequences')
    parser.add_argument('--show-products',
            help = 'pattern match and show show products.  For all products use empty string ("") (NOT IMPLEMENTED)')
    parser.add_argument('--region',
            metavar = 'start:end',
            type = lambda r: map(int, r.split(':')),
            help = 'parse specific region from genbank record')

def gb2info(seqname, seq, record):
    return {'seqname':seqname, 'tax_id':tax_of_genbank(record),
            'accession':record.id, 'description':record.description,
            'length':len(seq), 'ambig_count':count_ambiguous(seq)}

def action(args):
    records = islice(SeqIO.parse(args.infile, 'genbank'), args.limit)

    if args.type_strains:
        records = ifilter(lambda r: is_type(r), records)

    if args.filter:
        records = ifilter(lambda r: not UNCLASSIFIED_REGEX.search(r.description), records)

    info = []

    if args.features:
        args.features = set(args.features)
        # Parse out product locations
        for r in records:
            for f in r.features:
                products = set(f.qualifiers.get('product', []))
                if products & args.features:
                    tag = f.qualifiers.get('locus_tag', ['unspecified'])[0]
                    name = '{}_{}'.format(r.name, tag)
                    start, end = f.location.start.position, f.location.end.position
                    seq = r.seq[start:end]
                    length = len(seq)

                    if (args.minus and f.location.strand == 1) or f.location.strand == -1:
                        seq = seq.reverse_complement()

                    ambig_count = count_ambiguous(seq)

                    if length < args.min_length:
                        log.warning('dropping seq {} because of length {}'.format(name, length))
                        log.debug('Record and Feature information for short seq:')
                        log.debug(r)
                        log.debug(f)
                    elif ambig_count > args.max_ambiguous:
                        log.warning('dropping seq {} because of {} ambiguous bases'.format(
                            name, ambig_count))
                    else:
                        args.out.write('>{} {} {}\n{}\n'.format(name, r.id, r.description, seq))
                        info.append(gb2info(name, seq, r))
    else:
        # if no product specified output entire seq
        for r in records:
            if args.region:
                start, end = args.region[0], args.region[1]
                seq = r.seq[start:end]
                name = '{}_{}_{}'.format(r.name, start, end)
            else:
                seq = r.seq
                name = r.name

            length = len(seq)

            if length < args.min_length:
                log.warning('dropping seq {} because of length {}'.format(name, length))
                log.debug('Record and Feature information for short seq:')
                log.debug(r)
                log.debug(f)
            else:
                src = next(ifilter(lambda f: f.type == 'source', r.features), None)
                if src and ((args.minus and src.location.strand == 1) or src.location.strand == -1):
                    seq = seq.reverse_complement()

                args.out.write('>{} {} {}\n{}\n'.format(name, r.id, r.description, seq))
                info.append(gb2info(name, seq, r))

    if args.info_out:
        args.info_out.writeheader()
        args.info_out.writerows(info)

