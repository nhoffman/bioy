"""
Test classification.
"""

import filecmp
import logging
import sys

from os import path

import pandas as pd

from bioy_pkg.scripts.main import main

from __init__ import TestBase, TestCaseSuppressOutput, datadir as datadir

log = logging.getLogger(__name__)

class TestSubcommand(TestBase, TestCaseSuppressOutput):

    def main(self, arguments):
        main(['classifier'] + [str(a) for a in arguments])

    datadir = path.join(datadir, 'classifier')

    blast = path.join(datadir, 'rdp.blast.bz2')
    taxonomy = path.join(datadir, 'taxonomy.csv.bz2')
    seq_info = path.join(datadir, 'seq_info.csv.bz2')
    copy_numbers = path.join(datadir, 'rrnDB_16S_copy_num.csv')
    weights = path.join(datadir, 'weights.csv')

    def test01(self):
        """
        Minimal inputs.
        """

        outdir = self.mkoutdir()
        outfile = path.join(outdir, 'classifications.csv')

        self.main([
            '--seq-info', self.seq_info,
            '--taxonomy', self.taxonomy,
            '--outfile', outfile,
            self.blast])

        result = pd.read_csv(outfile)
        blast_data = pd.read_csv(self.blast, compression='bz2', header=None)
        self.assertEqual(result['reads'].sum(), len(blast_data[0].unique()))

    def test02(self):
        """
        Include weights.
        """

        outdir = self.mkoutdir()
        outfile = path.join(outdir, 'classifications.csv')

        self.main([
            '--seq-info', self.seq_info,
            '--taxonomy', self.taxonomy,
            '--weights', self.weights,
            '--outfile', outfile,
            self.blast])

        result = pd.read_csv(outfile)
        weights = pd.read_csv(self.weights, header=None)
        self.assertEqual(result['reads'].sum(), weights[1].sum())

