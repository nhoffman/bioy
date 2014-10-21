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

    def test01(self):
        """
        Test average basic usage
        """

        datadir = self.datadir

        this_test = sys._getframe().f_code.co_name
        sseqids = path.join(datadir, 'sseqids')

        outdir = self.mkoutdir()

        fasta_out = path.join(outdir, 'multi.fasta')
        info_out = path.join(outdir, 'seqinfo.csv')

        fasta_ref = path.join(datadir, 'multi.fasta')
        info_ref = path.join(datadir, 'seqinfo.csv')

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
        Ensure length of seqinfo matches fasta
        """

        datadir = self.datadir

        this_test = sys._getframe().f_code.co_name
        sseqids = path.join(datadir, 'sseqids')

        outdir = self.mkoutdir()

        fasta_out = path.join(outdir, 'multi.fasta')
        info_out = path.join(outdir, 'seqinfo.csv')

        fasta_ref = path.join(datadir, 'multi.fasta')
        info_ref = path.join(datadir, 'seqinfo.csv')

        args = ['--outfasta', fasta_out,
                '--seqinfo', info_out,
                '--email', 'ngh2@uw.edu',
                sseqids]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        fasta_lines = len([line for line in open(fasta_out)])
        info_lines = len([line for line in open(info_out)])

        self.assertTrue((info_lines-1)*2 == fasta_lines)
