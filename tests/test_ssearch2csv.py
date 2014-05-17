"""
Test ssearch2csv
"""

import filecmp
import logging

from os import path

from bioy_pkg.scripts.main import main

from __init__ import TestBase, TestCaseSuppressOutput, datadir

log = logging.getLogger(__name__)

class TestSsearch2csv(TestBase, TestCaseSuppressOutput):

    def main(self, arguments):
        main(['ssearch2csv'] + arguments)

    def test01(self):
        """
        Basic comparison of output
        """

        ssearch = path.join(datadir, 'rle_100_left.ssearch.bz2')
        outdir = self.mkoutdir()
        csv = path.join(outdir, 'ssearch.csv.bz2')
        self.main([ssearch, '--out', csv])

        reference = path.join(datadir, 'rle_100_left_ssearch.csv.bz2')

        self.assertTrue(filecmp.cmp(csv, reference))


