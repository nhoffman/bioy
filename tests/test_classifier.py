"""
Test classifier
"""

import logging

from os import path

import filecmp
import pandas as pd
import sys

from bioy_pkg.scripts.main import main

from __init__ import TestBase, TestCaseSuppressOutput, datadir as datadir

log = logging.getLogger(__name__)


class TestClassifier(TestBase, TestCaseSuppressOutput):

    def main(self, arguments):
        main(['classifier'] + [str(a) for a in arguments])

    log_info = 'bioy classify {}'

    datadir = path.join(datadir, 'classifier')
    thisdatadir = path.join(datadir, 'TestClassifier')

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

    def test03(self):
        """
        Test validation of type strains
        """

        thisdatadir = self.thisdatadir

        this_test = sys._getframe().f_code.co_name

        blast = path.join(thisdatadir, this_test, 'blast.csv.bz2')
        taxonomy = path.join(thisdatadir, this_test, 'taxonomy.csv.bz2')
        seq_info = path.join(thisdatadir, this_test, 'seq_info.csv.bz2')
        specimen_map = path.join(thisdatadir, this_test, 'map.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = path.join(outdir, 'classifications.csv.bz2')
        details_out = path.join(outdir, 'details.csv.bz2')

        classify_ref = path.join(
            thisdatadir, this_test, 'classifications.csv.bz2')
        details_ref = path.join(
            thisdatadir, this_test, 'details.csv.bz2')

        args = ['--specimen-map', specimen_map,
                '--details-out', details_out,
                '--out', classify_out,
                blast, seq_info, taxonomy]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))
