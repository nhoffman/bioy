"""
Test subcommands.
"""

from os import path
import logging
import os
import sys

from bioy_pkg.sequtils import fastalite
from bioy_pkg.scripts.main import main

from __init__ import TestBase, TestCaseSuppressOutput, datadir

log = logging.getLogger(__name__)

class TestDenoise(TestBase, TestCaseSuppressOutput):

    def main(self, arguments):
        """Pass list 'arguments' to scripts.main, ensuring that all are
        strings.

        """
        main(['denoise'] + [str(a) for a in arguments])


    def test01(self):
        fa = path.join(datadir, 'F1_3', 'trimmed.fasta')
        uc = path.join(datadir, 'F1_3', 'trimmed.uc')
        outdir = self.mkoutdir()
        fa_out = path.join(outdir, 'denoised.fasta')
        limit = 100
        self.main([fa, '--clusters', uc, '--outfile', fa_out, '--limit', limit])

        # cluster mass equals number of input sequences
        with open(fa_out) as f:
            cluster_mass = sum(int(line.split('_')[-1])
                               for line in f if line.startswith('>'))
            self.assertEqual(limit, cluster_mass)

        # regression test
        reference = path.join(datadir, 'F1_3', 'test01_denoised.fasta')
        with open(fa_out) as out, open(reference) as ref:
            outseqs = list(fastalite(out))
            refseqs = list(fastalite(ref))
            self.assertEqual(len(outseqs), len(refseqs))
            self.assertEqual(set(s.seq for s in outseqs), set(s.seq for s in refseqs))


    def test02(self):
        fa = path.join(datadir, 'F1_3', 'trimmed.fasta')
        uc = path.join(datadir, 'F1_3', 'trimmed.uc')
        outdir = self.mkoutdir()
        fa_out = path.join(outdir, 'denoised.fasta')
        limit = 100
        min_size = 2
        self.main([fa, '--clusters', uc, '--outfile', fa_out, '--limit', limit, '--min-clust-size',
                   min_size])

        reference = path.join(datadir, 'F1_3', 'test02_denoised.fasta')
        with open(fa_out) as out, open(reference) as ref:
            outseqs = list(fastalite(out))
            refseqs = list(fastalite(ref))
            self.assertEqual(len(outseqs), len(refseqs))
            self.assertEqual(set(s.seq for s in outseqs), set(s.seq for s in refseqs))


    def test03(self):
        fa = path.join(datadir, 'F1_3xR1_36', 'trimmed.fasta')
        uc = path.join(datadir, 'F1_3xR1_36', 'trimmed.uc')
        groups = path.join(datadir, 'F1_3xR1_36', 'groups.csv.bz2')
        outdir = self.mkoutdir()
        limit = 500
        min_size = 2

        denoised = path.join(outdir, 'denoised.fasta')
        self.main([fa,
                   '--clusters', uc,
                   '--outfile', denoised,
                   '--limit', limit,
                   '--min-clust-size', min_size,
                   '--weights', path.join(outdir, 'weights.csv')
               ])

        denoised_grouped = path.join(outdir, 'denoised.grouped.fasta')
        self.main([fa,
                   '--clusters', uc,
                   '--outfile', denoised_grouped,
                   '--limit', limit,
                   '--min-clust-size', min_size,
                   '--groups', groups,
                   '--weights', path.join(outdir, 'weights.grouped.csv')
               ])

        with open(denoised) as d, open(denoised_grouped) as g:
            ds = list(fastalite(d))
            gs = list(fastalite(g))
            self.assertEqual(set(s.seq for s in ds), set(s.seq for s in gs))

    def test04(self):
        """
        test no clusters passing min_size

        1) test file is actully created
        2) test there are no seqs in file
        """

        fa = path.join(datadir, 'F1_3', 'trimmed.fasta')
        uc = path.join(datadir, 'F1_3', 'trimmed.uc')
        outdir = self.mkoutdir()
        fa_out = path.join(outdir, 'denoised_empty.fasta')
        limit = 100
        min_size = sys.maxint
        self.main([fa, '--clusters', uc, '--outfile', fa_out, '--limit', limit, '--min-clust-size',
                   min_size])

        # 1)
        self.assertTrue(os.path.isfile(fa_out))

        with open(fa_out) as out:
            outseqs = list(fastalite(out))
            # 2)
            self.assertEqual(len(outseqs), 0)

    def test05(self):
        """
        test if no cluster file all one cluster
        """

        fa = path.join(datadir, '16S_random.fasta')
        outdir = self.mkoutdir()
        fa_out = fa_out = path.join(outdir, 'denoised.fasta')
        self.main([fa, '--outfile', fa_out])

        reference = path.join(datadir, '16S_random_cons.fasta')
        with open(fa_out) as out, open(reference) as ref:
            outseqs = list(fastalite(out))
            refseqs = list(fastalite(ref))
            self.assertEqual(len(outseqs), len(refseqs))
            self.assertEqual(set(s.seq for s in outseqs), set(s.seq for s in refseqs))

