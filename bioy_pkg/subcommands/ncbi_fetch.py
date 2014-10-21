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

"""
Fetch sequences from NCBI's nucleotide database using sequence identifiers (gi or gb)
Output is a multi-fasta of retrieved sequences and a corresponding sequence info (csv) file
"""

import logging
import sys

from Bio import Entrez
from Bio import SeqIO

from bioy_pkg.sequtils import FETCH_HEADERS
from bioy_pkg.utils import opener, Opener

def build_parser(parser):
    parser.add_argument('sseqids', nargs='?', type=Opener('r'),
            default = sys.stdin,
            help = 'input file, one identifier (gi or gb number) per line')
    parser.add_argument('-o', '--outfasta',
            type = Opener('w'),
            default = sys.stdout,
            help = 'multi-fasta, one sequence for each provided identifier')
    parser.add_argument('-i', '--seqinfo',
            type = Opener('w'),
            help = "optionally output seqinfo for each sequence : {}".format(FETCH_HEADERS))
    parser.add_argument('-h', '--no-header',
            help = "suppress seqinfo header")
    parser.add_argument('-e', '--email', required=True,
            help = "users of NCBI Entrez API should provide email.  if usage is excessive, ncbi may block access to its API")

def action(args):

    Entrez.email = args.email

    # cat bioyblast.out | cut -d, -f2 | cut -d\| -f2 > sseqids
    ids = [gi.strip() for gi in args.sseqids]

    handle = Entrez.efetch(db="nucleotide", id=','.join(ids), rettype="fasta", retmode="xml")
    records = Entrez.read(handle)
    
    if args.seqinfo and not args.no_header:
        args.seqinfo.write(','.join(FETCH_HEADERS) + '\n')

    # Note: some sequences may actually be proteins and use a larger alphabet than just ACTG
    for record in records:
        gi = record['TSeq_gi']
        gb = record['TSeq_accver']
        taxid = record['TSeq_taxid']
        orgname = record['TSeq_orgname']
        fullname = record['TSeq_defline']
        seq = record['TSeq_sequence']
        faheader = '|'.join(['>gi', gi, 'gb', gb, ''])+' '+fullname
        args.outfasta.write(faheader + '\n')
        args.outfasta.write(seq + '\n')
        if args.seqinfo:
            row = ','.join([gi, gb, taxid, orgname]) + '\n'
            args.seqinfo.write(row)
