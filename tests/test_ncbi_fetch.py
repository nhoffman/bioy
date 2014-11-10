"""
Test subcommands
"""

import filecmp
import logging
import sys

from os import path
from bioy_pkg.scripts.main import main

from __init__ import TestBase, TestCaseSuppressOutput, datadir as datadir

log = logging.getLogger(__name__)

class TestNcbiFetch(TestBase, TestCaseSuppressOutput):
    def main(self, arguments):
        main(['ncbi_fetch'] + [str(a) for a in arguments])

    log_info = 'bioy ncbi_fetch {}'

    datadir = path.join(datadir, 'ncbi_fetch')

###
# 01 - basic; whole sequence, no crop
# 02 - sequences cropped: from the middle, starting too early, ending too late, both
# 03 - non-existing sequence
###
    def test01(self):
        """
        Test basic usage
        """

        datadir = path.join(self.datadir)

        this_test = sys._getframe().f_code.co_name
        sseqids = path.join(datadir, 'sseqid_basic')

        outdir = self.mkoutdir()

        fasta_out = path.join(outdir, 'out1.fasta')
        info_out = path.join(outdir, 'out1.csv')

        fasta_ref = path.join(datadir, 'test01', 'ref.fasta')
        info_ref = path.join(datadir, 'test01', 'ref.csv')

        args = ['--outfasta', fasta_out,
                '--seqinfo', info_out,
                '--email', 'ngh2@uw.edu',
                sseqids]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(fasta_ref, fasta_out))
        self.assertTrue(filecmp.cmp(info_ref, info_out))

    def test02(self):
        """
        Cropped sequences
        """

        datadir = path.join(self.datadir)

        this_test = sys._getframe().f_code.co_name
        sseqids = path.join(datadir, 'sseqids_crops')

        outdir = self.mkoutdir()

        fasta_out = path.join(outdir, 'out.fasta')

        fasta_ref = path.join(datadir, 'test02', 'ref.fasta')

        args = ['--outfasta', fasta_out,
                '--email', 'ngh2@uw.edu',
                sseqids]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(fasta_ref, fasta_out))

    def test03(self):
        """
        ID with no matching sequence
        """

        datadir = path.join(self.datadir)

        this_test = sys._getframe().f_code.co_name
        sseqids = path.join(datadir, 'sseqid_missing')

        outdir = self.mkoutdir()

        fasta_out = path.join(outdir, 'out.fasta')

        args = ['--outfasta', fasta_out,
                '--email', 'ngh2@uw.edu',
                sseqids]

        log.info(self.log_info.format(' '.join(map(str, args))))

        # No output file is created.  Successfully passing the test means just 
        # gracefully handling the error without throwing an error

        self.main(args)
