"""
Test ssearch2csv
"""

import filecmp
import logging

from os import path

from bioy_pkg.subcommands import ssearch2csv

from __init__ import TestBase, TestSubcommand, datadir

log = logging.getLogger(__name__)

class TestSsearch2csv(TestBase, TestSubcommand):

    subcommand = ssearch2csv

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


