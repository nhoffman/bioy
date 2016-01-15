"""
Test subcommands
"""

import filecmp
import logging
import sys

from os import path
from bioy_pkg import main

from __init__ import TestBase, TestCaseSuppressOutput, datadir as datadir

log = logging.getLogger(__name__)

class TestNcbiFetch(TestBase, TestCaseSuppressOutput):
    def main(self, arguments):
        main(['ncbi_fetch'] + [str(a) for a in arguments])

    log_info = 'bioy ncbi_fetch {}'

    datadir = path.join(datadir, 'ncbi_fetch')

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
