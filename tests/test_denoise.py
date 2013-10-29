"""
Test subcommands.
"""

from os import path
import logging
import argparse

from bioy_pkg.utils import opener
from bioy_pkg.scripts.main import parse_arguments, main
from bioy_pkg.subcommands import denoise

from __init__ import TestCaseSuppressOutput, TestBase, TestSubcommand, datadir

log = logging.getLogger(__name__)


class TestDenoise(TestBase, TestSubcommand):

    subcommand = denoise

    def test01(self):
        fa = path.join(datadir, 'F1_3', 'trimmed.fasta')
        uc = path.join(datadir, 'F1_3', 'trimmed.uc')
        outdir = self.mkoutdir()
        fa_out = path.join(outdir, 'denoised.fasta')
        limit = 100
        self.main([fa, uc, '--outfile', fa_out, '--limit', limit])

        # cluster mass equals number of input sequences
        with open(fa_out) as f:
            cluster_mass = sum(int(line.split('_')[-1])
                               for line in f if line.startswith('>'))
            self.assertEqual(limit, cluster_mass)


    def test02(self):
        fa = path.join(datadir, 'F1_3', 'trimmed.fasta')
        uc = path.join(datadir, 'F1_3', 'trimmed.uc')
        outdir = self.mkoutdir()
        fa_out = path.join(outdir, 'denoised.fasta')
        limit = 100
        min_size = 2
        self.main([fa, uc, '--outfile', fa_out, '--limit', limit, '--min-clust-size',
                   min_size])

        # a regression test of sorts...
        reference = path.join(datadir, 'F1_3', 'test02_denoised.fasta')
        with open(fa_out) as out, open(reference) as ref:
            seqnames = lambda f: [line.strip()[1:] for line in f if line.startswith('>')]
            self.assertEqual(seqnames(out), seqnames(ref))
