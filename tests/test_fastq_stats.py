"""
Test fastq_stats subcommand
"""
import filecmp
import logging
from os import path

from bioy_pkg.scripts.main import main
from __init__ import TestBase, TestCaseSuppressOutput, datadir

log = logging.getLogger(__name__)

class TestFastq_stats(TestBase, TestCaseSuppressOutput):

    def main(self, arguments):
        main(['fastq_stats'] + arguments)

    fq_in = path.join(datadir, '16S_random.fastq')

    def test01(self):
        """
        Regression test for a simple test file.
        """

        outdir = self.mkoutdir()
        fq_out = path.join(outdir, 'fastq_stats.csv')
        reference = path.join(datadir, '16S_random.csv')
        args = ['--out', fq_out, self.fq_in]
        self.main(args)
        self.assertTrue(filecmp.cmp(reference, fq_out))

        #this_test = sys._getframe().f_code.co_name


