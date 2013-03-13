"""
Test cmscores subcommand
"""

from os import path
import unittest
import logging

from bioy_pkg.subcommands import cmscores

from . import TestBase, datadir

log = logging.getLogger(__name__)

class TestReadScores(TestBase):

    def test01(self):
        with open(path.join(datadir,'test1.cmscores')) as f:
            rows = list(cmscores.read_scores(f))
            self.assertEqual(len(rows[0]), len(cmscores.CMALIGN_SCORE_FIELDS))
            self.assertEqual(rows[0][0], '2')
            self.assertEqual(rows[-1][0], '195')

    def test02(self):
        with open(path.join(datadir,'test2.cmscores')) as f:
            rows = list(cmscores.read_scores(f))
            self.assertEqual(len(rows[0]), len(cmscores.CMALIGN_SCORE_FIELDS))
            self.assertEqual(rows[0][0], 'FUM0LCO01B97N3')
            self.assertEqual(rows[-1][0], 'FUM0LCO01D7IKL')
