"""
Test classifier
"""

import logging

from os import path

import filecmp
import sys

from bioy_pkg.scripts.main import main

from __init__ import TestBase, TestCaseSuppressOutput, datadir as datadir

log = logging.getLogger(__name__)


class TestClassifier(TestBase, TestCaseSuppressOutput):

    def main(self, arguments):
        main(['classifier'] + [str(a) for a in arguments])

    log_info = 'bioy classifier {}'

    copy_numbers = path.join(datadir, 'rrnDB_16S_copy_num.csv.bz2')

    thisdatadir = path.join(datadir, 'classifier', 'TestClassifier')

    def test01(self):
        """
        Minimal inputs.
        """

        this_test = sys._getframe().f_code.co_name

        thisdatadir = self.thisdatadir

        taxonomy = path.join(thisdatadir, 'taxonomy.csv.bz2')
        seq_info = path.join(thisdatadir, 'seq_info.csv.bz2')
        blast = path.join(thisdatadir, 'blast.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = path.join(outdir, 'classifications.csv.bz2')
        details_out = path.join(outdir, 'details.csv.bz2')

        classify_ref = path.join(
            thisdatadir, this_test, 'classifications.csv.bz2')
        details_ref = path.join(
            thisdatadir, this_test, 'details.csv.bz2')

        args = [
            '--out', classify_out,
            '--details-out', details_out,
            blast,
            seq_info,
            taxonomy]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))

    def test02(self):
        """
        Include weights.
        """

        this_test = sys._getframe().f_code.co_name

        thisdatadir = self.thisdatadir

        weights = path.join(thisdatadir, 'weights.csv.bz2')
        taxonomy = path.join(thisdatadir, 'taxonomy.csv.bz2')
        seq_info = path.join(thisdatadir, 'seq_info.csv.bz2')
        blast = path.join(thisdatadir, 'blast.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = path.join(outdir, 'classifications.csv.bz2')
        details_out = path.join(outdir, 'details.csv.bz2')

        classify_ref = path.join(
            thisdatadir, this_test, 'classifications.csv.bz2')
        details_ref = path.join(
            thisdatadir, this_test, 'details.csv.bz2')

        args = [
            '--weights', weights,
            '--out', classify_out,
            '--details-out', details_out,
            blast,
            seq_info,
            taxonomy]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))

    def test03(self):
        """
        Include specimen-map.
        """

        this_test = sys._getframe().f_code.co_name

        thisdatadir = self.thisdatadir

        specimen_map = path.join(thisdatadir, 'map.csv.bz2')
        taxonomy = path.join(thisdatadir, 'taxonomy.csv.bz2')
        seq_info = path.join(thisdatadir, 'seq_info.csv.bz2')
        blast = path.join(thisdatadir, 'blast.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = path.join(outdir, 'classifications.csv.bz2')
        details_out = path.join(outdir, 'details.csv.bz2')

        classify_ref = path.join(
            thisdatadir, this_test, 'classifications.csv.bz2')
        details_ref = path.join(
            thisdatadir, this_test, 'details.csv.bz2')

        args = [
            '--specimen-map', specimen_map,
            '--out', classify_out,
            '--details-out', details_out,
            blast,
            seq_info,
            taxonomy]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))

    def test05(self):
        """
        min-identity 99, max-identity 100
        """

        this_test = sys._getframe().f_code.co_name

        thisdatadir = self.thisdatadir

        taxonomy = path.join(thisdatadir, 'taxonomy.csv.bz2')
        seq_info = path.join(thisdatadir, 'seq_info.csv.bz2')
        blast = path.join(thisdatadir, 'blast.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = path.join(outdir, 'classifications.csv.bz2')
        details_out = path.join(outdir, 'details.csv.bz2')

        classify_ref = path.join(
            thisdatadir, this_test, 'classifications.csv.bz2')
        details_ref = path.join(
            thisdatadir, this_test, 'details.csv.bz2')

        args = [
            '--max-identity', '100',
            '--min-identity', '99',
            '--out', classify_out,
            '--details-out', details_out,
            blast,
            seq_info,
            taxonomy]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))

    def test06(self):
        """
        All together
        """

        this_test = sys._getframe().f_code.co_name

        thisdatadir = self.thisdatadir

        weights = path.join(thisdatadir, 'weights.csv.bz2')
        specimen_map = path.join(thisdatadir, 'map.csv.bz2')
        taxonomy = path.join(thisdatadir, 'taxonomy.csv.bz2')
        seq_info = path.join(thisdatadir, 'seq_info.csv.bz2')
        blast = path.join(thisdatadir, 'blast.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = path.join(outdir, 'classifications.csv.bz2')
        details_out = path.join(outdir, 'details.csv.bz2')

        classify_ref = path.join(
            thisdatadir, this_test, 'classifications.csv.bz2')
        details_ref = path.join(
            thisdatadir, this_test, 'details.csv.bz2')

        args = [
            '--max-identity', '100',
            '--min-identity', '99',
            '--specimen-map', specimen_map,
            '--weights', weights,
            '--copy-numbers', self.copy_numbers,
            '--out', classify_out,
            '--details-out', details_out,
            blast,
            seq_info,
            taxonomy]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))

    def test07(self):
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
    def test08(self):
        """
        Parse non-default blast result files that have haders
        """
        thisdatadir = self.thisdatadir

        this_test = sys._getframe().f_code.co_name

        blast = path.join(thisdatadir, 'blast_extrafields.csv.bz2')
        taxonomy = path.join(thisdatadir, 'taxonomy.csv.bz2')
        seq_info = path.join(thisdatadir, 'seq_info.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = path.join(outdir, 'classifications.csv.bz2')
        details_out = path.join(outdir, 'details.csv.bz2')

        classify_ref = path.join(
            thisdatadir, this_test, 'classifications.csv.bz2')
        details_ref = path.join(
            thisdatadir, this_test, 'details.csv.bz2')

        args = [
            '--has-header',
            '--out', classify_out,
            '--details-out', details_out,
            blast,
            seq_info,
            taxonomy]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))
