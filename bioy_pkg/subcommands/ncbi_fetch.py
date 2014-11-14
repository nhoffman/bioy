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

Note: If indexed bases are backwards (e.g. seq_stop > seq_stop) then they will be un-reversed
"""

import logging
import sys
import csv

from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import NucleotideAlphabet

from bioy_pkg.sequtils import FETCH_HEADERS
from bioy_pkg.utils import opener, Opener

fieldnames = ['id','seq_start','seq_stop']

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('sseqids', nargs='?', type=Opener('r'),
            default = sys.stdin,
            help = 'csv input file, each line containing (gi-or-gb),seq_start,seq_stop')
    parser.add_argument('-o', '--outfasta',
            type = Opener('w'),
            default = sys.stdout,
            help = 'multi-fasta, one sequence for each provided identifier')
    parser.add_argument('-i', '--seqinfo',
            type = Opener('w'),
            help = "optionally output seqinfo for each sequence : {}".format(FETCH_HEADERS))
    parser.add_argument('-n', '--no-header',
            help = "suppress seqinfo header")
    parser.add_argument('-e', '--email', required=True,
            help = "users of NCBI Entrez API should provide email.  if usage is excessive, ncbi may block access to its API")

def action(args):

    Entrez.email = args.email

    sseqids = csv.DictReader(args.sseqids, fieldnames=fieldnames)

    if args.seqinfo:
        info = csv.writer(args.seqinfo, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        if not args.no_header:
            info.writerow(FETCH_HEADERS)

    # For each subject line in the input, fetch the sequence and output to fasta
    for subject in sseqids:
        seq_handle = Entrez.efetch(db="nucleotide", id=subject['id'], 
                                   seq_start=subject['seq_start'], 
                                   seq_stop=subject['seq_stop'], 
                                   retmode="xml")

        sum_handle = Entrez.esummary(db="nucleotide", id=subject['id'])

        try:
            seq_records = Entrez.read(seq_handle)
            sum_records = Entrez.read(sum_handle)
        except:
            log.info("unable to parse seqid '{}'".format(subject['id']))
            continue

        assert(len(seq_records) == 1 and len(sum_records) == 1)
            
        seq_record = seq_records[0]
        sum_record = sum_records[0]

        # Ideally we would assert that the seqids of the two match, but one uses gi and the other uses gb

        taxid = sum_record['TaxId']
        description = seq_record['GBSeq_definition']
        sequence = seq_record['GBSeq_sequence'].upper()
        if subject['seq_start'] and subject['seq_stop']:
            assert(len(sequence) <= abs(int(subject['seq_stop']) - int(subject['seq_start']))+1)

        seqids = [seqid for seqid in seq_record['GBSeq_other-seqids'] if 'gi' in seqid or 'gb' in seqid]
        assert(len(seqids) >= 1)
        seqid = ''.join(seqids) # Joins genbank, gi numbers, etc, which include '|'
        outseq = SeqRecord(Seq(sequence, NucleotideAlphabet),
                           id=seqid,
                           description=description)
        SeqIO.write(outseq, args.outfasta, 'fasta')
        if args.seqinfo:
            info.writerow([seqid, taxid, description])
    args.outfasta.close()
