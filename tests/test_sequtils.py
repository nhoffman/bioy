"""
Test sequtils module.
"""

from os import path
import unittest
import logging
from itertools import chain

from bioy_pkg import sequtils

from __init__ import TestBase
log = logging.getLogger(__name__)

class TestRunSsearch(TestBase):

    def test01(self):
        """
        Provide an output name
        """

        q,t = self.data('two.fasta'), self.data('ten.fasta')
        out = path.join(self.mkoutdir(), 'aligns.ssearch')
        with sequtils.run_ssearch(q, t, out) as aligns:
            self.assertEqual(out, aligns.name)
            parsed = sequtils.parse_ssearch36(aligns)
            self.assertEqual(set(['H59735', 'T70875']),
                             {d['q_name'] for d in parsed})

        # should still exist since we provided a name for the output
        self.assertTrue(path.exists(out))

    def test02(self):
        """
        No output name, don't clean up
        """

        q,t = self.data('two.fasta'), self.data('ten.fasta')
        with sequtils.run_ssearch(q, t, cleanup=False) as aligns:
            parsed = sequtils.parse_ssearch36(aligns)
            self.assertEqual(set(['H59735', 'T70875']),
                             {d['q_name'] for d in parsed})

        # should still exist since we specified 'cleanup=False'
        self.assertTrue(path.exists(aligns.name))

    def test03(self):
        """
        No output name, clean up
        """

        q,t = self.data('two.fasta'), self.data('ten.fasta')
        with sequtils.run_ssearch(q, t) as aligns:
            parsed = sequtils.parse_ssearch36(aligns)
            self.assertEqual(set(['H59735', 'T70875']),
                             {d['q_name'] for d in parsed})

        # should not still exist
        self.assertFalse(path.exists(aligns.name))

class TestAllPairwise(TestBase):

    def test01(self):
        with open(self.data('five.fasta')) as f:
            seqs = list(sequtils.fastalite(f))
            pairs = list(sequtils.all_pairwise(seqs))
            self.assertEqual(len(pairs), (len(seqs)*(len(seqs)-1))/2)
            self.assertEqual([s.id for s in seqs], list(sequtils.names_from_pairs(pairs)))

    def test02(self):
        with open(self.data('two.fasta')) as f:
            seqs = list(sequtils.fastalite(f))
            pairs = list(sequtils.all_pairwise(seqs))
            self.assertEqual(len(pairs), (len(seqs)*(len(seqs)-1))/2)
            self.assertEqual([s.id for s in seqs], list(sequtils.names_from_pairs(pairs)))