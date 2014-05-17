"""
Test subcommands.
"""

import filecmp
import logging

from os import path

from bioy_pkg.scripts.main import main

from __init__ import TestBase, TestCaseSuppressOutput, datadir as datadir

log = logging.getLogger(__name__)

class TestCsv2fasta(TestBase, TestCaseSuppressOutput):

    log_info = 'bioy csv2fasta {}'

    datadir = path.join(datadir, 'csv2fasta')

    def main(self, arguments):
        main(['csv2fasta'] + arguments)

    def test01(self):
        """
        Test construction of fasta file from ssearch results
        """

        datadir = self.datadir

        outdir = self.mkoutdir()

        ssearch_in = path.join(datadir, 'ssearch.csv.bz2')

        fa_out = path.join(outdir, 'ssearch.fasta.bz2')

        reference = path.join(datadir, 'ssearch.fasta.bz2')

        args = ['--columns', 'q_name,q_seq', '--out', fa_out, ssearch_in]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(reference, fa_out))

