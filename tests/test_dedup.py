"""
Test dedup subcommand
"""

import filecmp
import logging
import sys

from os import path

from bioy_pkg import main

from __init__ import TestBase, TestCaseSuppressOutput, datadir as datadir

log = logging.getLogger(__name__)

class TestDedup(TestBase, TestCaseSuppressOutput):

    def main(self, arguments):
        main(['dedup'] + arguments)

    log_info = 'bioy dedup {}'

    datadir = path.join(datadir, 'dedup')

    fa_in = path.join(datadir, 'seqs.fasta.bz2')

    split_info = path.join(datadir, 'split_info.csv.bz2')

    def test01(self):
        """
        Test basic usage with no split_info
        """

        datadir = self.datadir

        outdir = self.mkoutdir()

        this_test = sys._getframe().f_code.co_name

        dedup_out = path.join(outdir, 'dedup.fasta.bz2')

        reference = path.join(datadir, this_test, 'dedup.fasta.bz2')

        args = ['--out', dedup_out, self.fa_in]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(reference, dedup_out))

    def test02(self):
        """
        Test primary and secondary group species and tax_id
        """

        datadir = self.datadir

        outdir = self.mkoutdir()

        this_test = sys._getframe().f_code.co_name

        dedup_out = path.join(outdir, 'dedup.fasta.bz2')

        reference = path.join(datadir, this_test, 'dedup.fasta.bz2')

        args = ['--primary-group', 'species', '--secondary-group', 'tax_id',
                '--split-info', self.split_info, '--out', dedup_out, self.fa_in]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(reference, dedup_out))

    def test03(self):
        """
        Test primary and secondary group phylum and tax_id
        """

        datadir = self.datadir

        outdir = self.mkoutdir()

        this_test = sys._getframe().f_code.co_name

        dedup_out = path.join(outdir, 'dedup.fasta.bz2')

        reference = path.join(datadir, this_test, 'dedup.fasta.bz2')

        args = ['--primary-group', 'phylum', '--secondary-group', 'tax_id',
                '--split-info', self.split_info, '--out', dedup_out, self.fa_in]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(reference, dedup_out))

    def test04(self):
        """
        Test weights out with primary and secondary group species and tax_id
        """

        datadir = self.datadir

        outdir = self.mkoutdir()

        this_test = sys._getframe().f_code.co_name

        dedup_out = path.join(outdir, 'dedup.fasta.bz2')
        weights_out = path.join(outdir, 'weights.fasta.bz2')

        dedup_ref = path.join(datadir, this_test, 'dedup.fasta.bz2')
        weights_ref = path.join(datadir, this_test, 'weights.fasta.bz2')

        args = ['--primary-group', 'species', '--secondary-group', 'tax_id',
                '--out-weights', weights_out,
                '--split-info', self.split_info, '--out', dedup_out,
                self.fa_in]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(dedup_ref, dedup_out))
        self.assertTrue(filecmp.cmp(weights_ref, weights_out))

    def test05(self):
        """
        Test weights out with no split-info and primary and secondary group species and tax_id
        """

        datadir = self.datadir

        outdir = self.mkoutdir()

        this_test = sys._getframe().f_code.co_name

        dedup_out = path.join(outdir, 'dedup.fasta.bz2')
        weights_out = path.join(outdir, 'weights.fasta.bz2')

        dedup_ref = path.join(datadir, this_test, 'dedup.fasta.bz2')
        weights_ref = path.join(datadir, this_test, 'weights.fasta.bz2')

        args = ['--primary-group', 'species', '--secondary-group', 'tax_id',
                '--out-weights', weights_out, '--out', dedup_out, self.fa_in]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(dedup_ref, dedup_out))
        self.assertTrue(filecmp.cmp(weights_ref, weights_out))

    def test06(self):
        """
        Test weights out with no split-info and primary and secondary group species and tax_id
        """

        datadir = self.datadir

        outdir = self.mkoutdir()

        this_test = sys._getframe().f_code.co_name

        dedup_out = path.join(outdir, 'dedup.fasta.bz2')
        map_out = path.join(outdir, 'map.fasta.bz2')

        dedup_ref = path.join(datadir, this_test, 'dedup.fasta.bz2')
        map_ref = path.join(datadir, this_test, 'map.fasta.bz2')

        args = ['--primary-group', 'species', '--secondary-group', 'tax_id',
                '--out-map', map_out, '--out', dedup_out, self.fa_in]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(dedup_ref, dedup_out))
        self.assertTrue(filecmp.cmp(map_ref, map_out))

    def test07(self):
        """
        Test everything together with primary and secondary group species and tax_id
        """

        datadir = self.datadir

        outdir = self.mkoutdir()

        this_test = sys._getframe().f_code.co_name

        dedup_out = path.join(outdir, 'dedup.fasta.bz2')
        map_out = path.join(outdir, 'map.fasta.bz2')
        weights_out = path.join(outdir, 'weights.fasta.bz2')

        dedup_ref = path.join(datadir, this_test, 'dedup.fasta.bz2')
        map_ref = path.join(datadir, this_test, 'map.fasta.bz2')
        weights_ref = path.join(datadir, this_test, 'weights.fasta.bz2')

        args = ['--primary-group', 'species', '--secondary-group', 'tax_id',
                '--split-info', self.split_info,
                '--out-weights', weights_out,
                '--out-map', map_out,
                '--out', dedup_out, self.fa_in]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(dedup_ref, dedup_out))
        self.assertTrue(filecmp.cmp(map_ref, map_out))
        self.assertTrue(filecmp.cmp(weights_ref, weights_out))

