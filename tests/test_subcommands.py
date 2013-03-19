"""
Test subcommands.
"""

import argparse
import logging

from bioy_pkg.scripts.main import main
from __init__ import TestCaseSuppressOutput, TestBase

log = logging.getLogger(__name__)

class TestMain(TestCaseSuppressOutput, TestBase):

    def testExit01(self):
        self.assertRaises(SystemExit, main, ['-h'])
