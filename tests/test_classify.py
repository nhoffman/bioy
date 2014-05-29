"""
Test subcommands.
"""

import filecmp
import logging
import sys

from os import path

from bioy_pkg.scripts.main import main

from __init__ import TestBase, TestCaseSuppressOutput, datadir as datadir

log = logging.getLogger(__name__)

class TestClassify(TestBase, TestCaseSuppressOutput):

    def main(self, arguments):
        main(['classify'] + [str(a) for a in arguments])

    log_info = 'bioy classify {}'

    copy_numbers = path.join(datadir, 'rrnDB_16S_copy_num.csv.bz2')

    datadir = path.join(datadir, 'classify')

    tax = path.join(datadir, 'taxonomy.csv.bz2')
    info = path.join(datadir, 'seq_info.csv.bz2')

    def test01(self):
        """
        Test average basic usage
        """

        datadir = self.datadir

        this_test = sys._getframe().f_code.co_name

        blast = path.join(datadir, 'blast.csv')
        weights = path.join(datadir, 'weights.csv.bz2')
        mapp = path.join(datadir, 'map.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = path.join(outdir, 'classification.csv.bz2')
        details_out = path.join(outdir, 'details.csv.bz2')

        classify_ref = path.join(datadir, this_test, 'classification.csv.bz2')
        details_ref = path.join(datadir, this_test, 'details.csv.bz2')

        args = ['--map', mapp,
                '--weights', weights,
                '--seq-info', self.info,
                '--tax', self.tax,
                '--copy-numbers', self.copy_numbers,
                '--max-identity', 100,
                '--min-identity', 90,
                '--coverage', 90,
                '--asterisk', 100,
                '--target-max-group', 3,
                '--group-def', 'species:4',
                '--target-rank', 'species',
                '--details-identity', 90,
                '--out-detail', details_out,
                '--out', classify_out,
                blast]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))

    def test02(self):
        """
        Test no mapp
        """

        datadir = self.datadir

        this_test = sys._getframe().f_code.co_name

        blast = path.join(datadir, 'blast.csv')
        weights = path.join(datadir, 'weights.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = path.join(outdir, 'classification.csv.bz2')
        details_out = path.join(outdir, 'details.csv.bz2')

        classify_ref = path.join(datadir, this_test, 'classification.csv.bz2')
        details_ref = path.join(datadir, this_test, 'details.csv.bz2')

        args = ['--weights', weights,
                '--seq-info', self.info,
                '--tax', self.tax,
                '--copy-numbers', self.copy_numbers,
                '--max-identity', 100,
                '--min-identity', 90,
                '--coverage', 90,
                '--asterisk', 100,
                '--target-max-group', 3,
                '--group-def', 'species:4',
                '--target-rank', 'species',
                '--details-identity', 90,
                '--out-detail', details_out,
                '--out', classify_out,
                blast]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))

    def test03(self):
        """
        Test no weights
        """

        datadir = self.datadir

        this_test = sys._getframe().f_code.co_name

        blast = path.join(datadir, 'blast.csv')
        mapp = path.join(datadir, 'map.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = path.join(outdir, 'classification.csv.bz2')
        details_out = path.join(outdir, 'details.csv.bz2')

        classify_ref = path.join(datadir, this_test, 'classification.csv.bz2')
        details_ref = path.join(datadir, this_test, 'details.csv.bz2')

        args = ['--map', mapp,
                '--seq-info', self.info,
                '--tax', self.tax,
                '--copy-numbers', self.copy_numbers,
                '--max-identity', 100,
                '--min-identity', 90,
                '--coverage', 90,
                '--asterisk', 100,
                '--target-max-group', 3,
                '--group-def', 'species:4',
                '--target-rank', 'species',
                '--details-identity', 90,
                '--out-detail', details_out,
                '--out', classify_out,
                blast]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))

    def test04(self):
        """
        Test no weights, no map
        """

        datadir = self.datadir

        this_test = sys._getframe().f_code.co_name

        blast = path.join(datadir, 'blast.csv')

        outdir = self.mkoutdir()

        classify_out = path.join(outdir, 'classification.csv.bz2')
        details_out = path.join(outdir, 'details.csv.bz2')

        classify_ref = path.join(datadir, this_test, 'classification.csv.bz2')
        details_ref = path.join(datadir, this_test, 'details.csv.bz2')

        args = ['--seq-info', self.info,
                '--tax', self.tax,
                '--copy-numbers', self.copy_numbers,
                '--max-identity', 100,
                '--min-identity', 90,
                '--coverage', 90,
                '--asterisk', 100,
                '--target-max-group', 3,
                '--group-def', 'species:4',
                '--target-rank', 'species',
                '--details-identity', 90,
                '--out-detail', details_out,
                '--out', classify_out,
                blast]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))

    def test05(self):
        """
        Test grouping by 1
        """

        datadir = self.datadir

        this_test = sys._getframe().f_code.co_name

        blast = path.join(datadir, 'blast.csv')
        weights = path.join(datadir, 'weights.csv.bz2')
        mapp = path.join(datadir, 'map.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = path.join(outdir, 'classification.csv.bz2')
        details_out = path.join(outdir, 'details.csv.bz2')

        classify_ref = path.join(datadir, this_test, 'classification.csv.bz2')
        details_ref = path.join(datadir, this_test, 'details.csv.bz2')

        args = ['--map', mapp,
                '--weights', weights,
                '--seq-info', self.info,
                '--tax', self.tax,
                '--copy-numbers', self.copy_numbers,
                '--max-identity', 100,
                '--min-identity', 90,
                '--coverage', 90,
                '--asterisk', 100,
                '--target-max-group', 1,
                '--group-def', 'species:1',
                '--target-rank', 'species',
                '--details-identity', 90,
                '--out-detail', details_out,
                '--out', classify_out,
                blast]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))

    def test06(self):
        """
        Test high level grouping
        """

        datadir = self.datadir

        this_test = sys._getframe().f_code.co_name

        blast = path.join(datadir, 'blast.csv')
        weights = path.join(datadir, 'weights.csv.bz2')
        mapp = path.join(datadir, 'map.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = path.join(outdir, 'classification.csv.bz2')
        details_out = path.join(outdir, 'details.csv.bz2')

        classify_ref = path.join(datadir, this_test, 'classification.csv.bz2')
        details_ref = path.join(datadir, this_test, 'details.csv.bz2')

        args = ['--map', mapp,
                '--weights', weights,
                '--seq-info', self.info,
                '--tax', self.tax,
                '--copy-numbers', self.copy_numbers,
                '--max-identity', 100,
                '--min-identity', 90,
                '--coverage', 90,
                '--asterisk', 100,
                '--target-max-group', 1,
                '--group-def', 'species:1',
                '--target-rank', 'superkingdom',
                '--details-identity', 90,
                '--out-detail', details_out,
                '--out', classify_out,
                blast]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))

    def test07(self):
        """
        Test not copy number correction
        """

        datadir = self.datadir

        this_test = sys._getframe().f_code.co_name

        blast = path.join(datadir, 'blast.csv')
        weights = path.join(datadir, 'weights.csv.bz2')
        mapp = path.join(datadir, 'map.csv.bz2')

        outdir = self.mkoutdir()

        classify_out = path.join(outdir, 'classification.csv.bz2')
        details_out = path.join(outdir, 'details.csv.bz2')

        classify_ref = path.join(datadir, this_test, 'classification.csv.bz2')
        details_ref = path.join(datadir, this_test, 'details.csv.bz2')

        args = ['--map', mapp,
                '--weights', weights,
                '--seq-info', self.info,
                '--tax', self.tax,
                '--max-identity', 100,
                '--min-identity', 90,
                '--coverage', 90,
                '--asterisk', 100,
                '--target-max-group', 3,
                '--group-def', 'species:4',
                '--target-rank', 'species',
                '--details-identity', 90,
                '--out-detail', details_out,
                '--out', classify_out,
                blast]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(classify_ref, classify_out))
        self.assertTrue(filecmp.cmp(details_ref, details_out))

