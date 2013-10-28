"""
Test subcommands.
"""

from os import path
import logging
import argparse

from bioy_pkg.utils import opener
from bioy_pkg.scripts.main import parse_arguments, main
from bioy_pkg.subcommands import reverse_complement

from __init__ import TestCaseSuppressOutput, TestBase, TestSubcommand, datadir

log = logging.getLogger(__name__)


def make_caller(subcommand):
    """Factory function to construct a test case method for calling the
    subcommand's action.

    """

    def _parser(_args):
        parser = argparse.ArgumentParser()
        subcommand.build_parser(parser)
        return subcommand.action(parser.parse_args(_args))

    return _parser

class TestMainScript(TestCaseSuppressOutput, TestBase):

    def testExit01(self):
        self.assertRaises(SystemExit, main, ['notacommand'])

    def testExit02(self):
        self.assertRaises(SystemExit, main, ['-h'])

class TestReverseComplement(TestBase, TestSubcommand):

    subcommand = reverse_complement

    def test01(self):
        outdir = self.mkoutdir()
        fa = path.join(datadir, 'F1_3', 'trimmed_rle.fasta')
        fa_out = path.join(outdir, 'rc.fasta')
        self.main([fa, '-o', fa_out])
        self.assertTrue(path.exists(fa_out))

        with open(fa) as infile, open(fa_out) as outfile:
            nseqs = lambda f: sum(1 for line in f if line.startswith('>'))
            self.assertEqual(nseqs(infile), nseqs(outfile))

    def test02(self):
        outdir = self.mkoutdir()
        fa = path.join(datadir, 'F1_3', 'trimmed_rle.fasta')
        rle = path.join(datadir, 'F1_3', 'trimmed_rle.csv.bz2')

        fa_out = path.join(outdir, 'rc.fasta')
        rle_out = path.join(outdir, 'rc.csv.bz2')

        self.main([fa, rle, '--out-fasta', fa_out, '--out-rle', rle_out])
        self.assertTrue(path.exists(fa_out))
        self.assertTrue(path.exists(rle_out))

        with opener(rle) as infile, opener(rle_out) as outfile:
            self.assertEqual(infile.next(), outfile.next())
            self.assertEqual(len(list(infile)), len(list(outfile)))
