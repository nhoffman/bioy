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
Fetch nucleotide sequences from NCBI using sequence identifiers
"""

import logging
import sys

#from csv import DictWriter
#from cStringIO import StringIO
#from itertools import chain, groupby
#from operator import itemgetter
#from subprocess import Popen, PIPE

from Bio import Entrez
from Bio import SeqIO

#from bioy_pkg.sequtils import BLAST_HEADER, BLAST_FORMAT, fastalite
from bioy_pkg.utils import opener, Opener

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('sseqids', nargs='?', type=Opener('r'),
            default = sys.stdin,
            help = 'input file, one gi per line')
    parser.add_argument('-o', '--outfasta',
            type = Opener('w'),
            default = sys.stdout,
            help = 'Multi-fasta')
    parser.add_argument('-i', '--seqinfo',
            type= Opener('w'),
            help = "Optionally output seqinfo for each sequence")

def action(args):

    Entrez.email = "ngh2@uw.edu"

    # cat bioyblast.out | cut -d, -f2 | cut -d\| -f2 > sseqids
    ids = [gi.strip() for gi in args.sseqids]

    handle = Entrez.efetch(db="nucleotide", id=','.join(ids), rettype="fasta", retmode="xml")
    records = Entrez.read(handle)
    
    if args.seqinfo:
        args.seqinfo.write(','.join(['gi','gb','taxid','orgname']) + '\n')
    for record in records:
        if record['TSeq_seqtype'].attributes['value'] != 'nucleotide':
            if record['TSeq_seqtype'].attributes['value'] not in ['protein']:
                print("Skipping {} record".format(record['TSeq_seqtype'].attributes['value']))
            continue
        gi = record['TSeq_gi']
        gb = record['TSeq_accver']
        taxid = record['TSeq_taxid']
        orgname = record['TSeq_orgname']
        fullname = record['TSeq_defline']
        seq = record['TSeq_sequence']
        faheader = '|'.join(['>gi', gi, 'gb', gb, ''])+' '+fullname
        args.outfasta.write(faheader+'\n')
        args.outfasta.write(seq+'\n')
        if args.seqinfo:
            row = ','.join([gi, gb, taxid, orgname]) + '\n'
            args.seqinfo.write(row)
