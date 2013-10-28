import sys
import logging
import os
from os import path
import unittest
import argparse

from bioy_pkg.utils import mkdir

# set up logging for unit tests
verbosity_flag = [x for x in sys.argv if x.startswith('-v')]
verbosity = (verbosity_flag[0] if verbosity_flag else '').count('v')

loglevel = {
    0: logging.WARNING,
    1: logging.INFO,
    2: logging.DEBUG,
}.get(verbosity, logging.DEBUG)

if verbosity > 1:
    logformat = '%(levelname)s %(module)s %(lineno)s %(message)s'
else:
    logformat = '%(message)s'

logging.basicConfig(file=sys.stdout, format=logformat, level=loglevel)
log = logging.getLogger(__name__)

# module data
datadir = 'testfiles'
outputdir = 'test_output'

mkdir(outputdir)

class TestBase(unittest.TestCase):
    """
    Base class for unit tests with methods for defining output
    directories based on method name.
    """

    outputdir = outputdir

    def mkoutdir(self, clobber = True):
        """
        Create outdir as outpudir/module.class.method (destructively
        if clobber is True).
        """

        funcname = '.'.join(self.id().split('.')[-3:])
        outdir = path.join(self.outputdir, funcname)
        mkdir(outdir, clobber)
        return outdir

    def data(self, fname):
        return path.join(datadir, fname)


class TestCaseSuppressOutput(unittest.TestCase):

    def setUp(self):
        self.funcname = '_'.join(self.id().split('.')[-2:])
        self.suppress_output = log.getEffectiveLevel() >= logging.INFO
        if self.suppress_output:
            sys.stdout = sys.stderr = open(os.devnull, 'w')

    def tearDown(self):
        if self.suppress_output:
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__

class TestSubcommand(unittest.TestCase):
    """Must define class variable subcommand with methods action and
    build_parser. Execute the subcommand using

      self.main(list_of_args)

    """

    def setUp(self):
        self.funcname = '_'.join(self.id().split('.')[-2:])
        self.suppress_output = log.getEffectiveLevel() >= logging.INFO
        if self.suppress_output:
            sys.stdout = sys.stderr = open(os.devnull, 'w')

        self.action = self.subcommand.action
        self.build_parser = self.subcommand.build_parser
        parser = argparse.ArgumentParser()
        self.build_parser(parser)
        self.main = lambda args: self.action(parser.parse_args(args))


    def tearDown(self):
        if self.suppress_output:
            sys.stdout = sys.__stdout__
            sys.stderr = sys.__stderr__
