"""
Test subcommands.
"""

import filecmp
import logging

from os import path

from bioy_pkg.scripts.main import main

from __init__ import TestBase, TestCaseSuppressOutput, datadir as datadir

log = logging.getLogger(__name__)

class TestSsearch_count(TestBase, TestCaseSuppressOutput):

    def main(self, arguments):
        main(['ssearch_count'] + [str(a) for a in arguments])

    log_info = 'bioy ssearch_count {}'

    datadir = path.join(datadir, 'ssearch_count')

    tax = path.join(datadir, 'taxonomy.csv.bz2')
    info = path.join(datadir, 'seq_info.csv.bz2')

    def test01(self):
        """
        Test average global superkingdom breakdown with zscore
        filter of 100
        """

        datadir = self.datadir

        ssearch = path.join(datadir, 'ssearch.csv.bz2')

        outdir = self.mkoutdir()

        out = path.join(outdir, 'superkingdom.csv.bz2')

        ref = path.join(datadir, 'test01', 'superkingdom.csv.bz2')

        args = ['--min-zscore', '100',
                '--taxonomy', self.tax,
                '--info', self.info,
                '--rank', 'superkingdom',
                '--out', out,
                ssearch]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(ref, out))

    def test02(self):
        """
        Test average global species breakdown with no zscore
        """

        datadir = self.datadir

        ssearch = path.join(datadir, 'ssearch.csv.bz2')

        outdir = self.mkoutdir()

        out = path.join(outdir, 'superkingdom.csv.bz2')

        ref = path.join(datadir, 'test02', 'superkingdom.csv.bz2')

        args = ['--taxonomy', self.tax,
                '--info', self.info,
                '--rank', 'superkingdom',
                '--out', out,
                ssearch]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(ref, out))

    def test03(self):
        """
        Test average species breakdown with zscore 100
        """

        datadir = self.datadir

        ssearch = path.join(datadir, 'ssearch.csv.bz2')

        outdir = self.mkoutdir()

        out = path.join(outdir, 'species.csv.bz2')

        ref = path.join(datadir, 'test03', 'species.csv.bz2')

        args = ['--min-zscore', '100',
                '--taxonomy', self.tax,
                '--info', self.info,
                '--rank', 'species',
                '--out', out,
                ssearch]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(ref, out))

    def test04(self):
        """
        Test average class breakdown with zscore 100
        """

        datadir = self.datadir

        ssearch = path.join(datadir, 'ssearch.csv.bz2')

        outdir = self.mkoutdir()

        out = path.join(outdir, 'class.csv.bz2')

        ref = path.join(datadir, 'test04', 'class.csv.bz2')

        args = ['--min-zscore', '100',
                '--taxonomy', self.tax,
                '--info', self.info,
                '--rank', 'class',
                '--out', out,
                ssearch]

        log.info(self.log_info.format(' '.join(map(str, args))))

        self.main(args)

        self.assertTrue(filecmp.cmp(ref, out))

