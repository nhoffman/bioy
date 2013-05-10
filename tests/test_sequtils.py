"""
Test sequtils module.
"""

import logging
import pprint

from os import path

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

class TestTempFasta(TestBase):

    def setUp(self):
        self.outdir = self.mkoutdir()
        with open(self.data('two.fasta')) as f:
            self.seqs = list(sequtils.fastalite(f))

    def test01(self):
        with sequtils.fasta_tempfile(self.seqs) as f:
            self.assertTrue(path.exists(f))
        self.assertFalse(path.exists(f))

    def test02(self):
        with sequtils.fasta_tempfile(self.seqs, dir = self.outdir) as f:
            self.assertTrue(path.exists(f))
        self.assertFalse(path.exists(f))

class TestEncodeAndDecode(TestBase):
    def test01(self):
        seq = 'TCTGGACCGTGTCTTTCAGTTCCAAAGTGTGACTGATCCATCCTCTCAGACC'

        e,c = sequtils.homoencode(seq)

        self.assertEquals(seq, sequtils.homodecode(e,c))

class TestRle(TestBase):

    def test01(self):
        rle,counts = sequtils.homoencode('AAA')
        self.assertEquals(rle,'A')
        self.assertEquals(counts,[3])

    def test02(self):
        rle,counts = sequtils.homoencode('GCTTCAAACATA')
        self.assertEquals(rle, 'GCTCACATA')
        self.assertEquals(counts, [1,1,2,1,3,1,1,1,1])

class TestDecodeAlignment(TestBase):

    def test01(self):
        pass

class TestErrorCounting(TestBase):

    def tearDown(self):

        # print '\nref:  ', self.ref
        # print 'query:', self.query
        errors = list(sequtils.itemize_errors(self.ref, self.query))
        log.debug(sequtils.show_errors(errors))

        log.debug('\n'+pprint.pformat(self.expected))
        log.debug('\n'+pprint.pformat(errors))

        # indices into ref are as expected
        self.assertEquals(set(e['i'] for e in errors), set(self.expected.keys()))

        # individual positions are as expected
        for pos in errors:
            self.assertEquals(pos['ref'], self.expected[pos['i']]['ref'])
            self.assertEquals(pos['query'], self.expected[pos['i']]['query'])

    def test02(self):

        #             01234-56789
        self.ref =   'GATTA-CATA-'
        self.query = 'GCTTAACATAG'

        self.expected = {
            0: dict(ref='G',query='G'),
            1: dict(ref='A',query='C'),
            2: dict(ref='TT',query='TT'),
            4: dict(ref='A-',query='AA'),
            5: dict(ref='C',query='C'),
            6: dict(ref='A',query='A'),
            7: dict(ref='T',query='T'),
            8: dict(ref='A',query='A')
            }

    def test03(self):

        #             01234-56789
        self.ref =   'GATT-ACATA-'
        self.query = 'GCTTAACATAG'

        self.expected = {
            0: dict(ref='G',query='G'),
            1: dict(ref='A',query='C'),
            2: dict(ref='TT',query='TT'),
            4: dict(ref='-A',query='AA'),
            5: dict(ref='C',query='C'),
            6: dict(ref='A',query='A'),
            7: dict(ref='T',query='T'),
            8: dict(ref='A',query='A')
            }

    def test04(self):

        #             01234567890
        self.ref =   'GCTTAACATAG'
        self.query = 'GATT-ACATA-'

        self.expected = {
            0: dict(ref='G',query='G'),
            1: dict(ref='C',query='A'),
            2: dict(ref='TT',query='TT'),
            4: dict(ref='AA',query='-A'),
            6: dict(ref='C',query='C'),
            7: dict(ref='A',query='A'),
            8: dict(ref='T',query='T'),
            9: dict(ref='A',query='A'),
            10: dict(ref='G',query='-')
            }

    def test05(self):

        #              0123456789
        self.ref =   '-CTTAACATAG'
        self.query = 'GATT-ACATAG'

        self.expected = {
            0: dict(ref='C',query='A'),
            1: dict(ref='TT',query='TT'),
            3: dict(ref='AA',query='-A'),
            5: dict(ref='C',query='C'),
            6: dict(ref='A',query='A'),
            7: dict(ref='T',query='T'),
            8: dict(ref='A',query='A'),
            9: dict(ref='G',query='G')
            }

    def test06(self):

        #             01234567890
        self.ref =   'GCTTAACATAG'
        self.query = '-ATTA=CATA-'

        self.expected = {
            0: dict(ref='G',query='-'),
            1: dict(ref='C',query='A'),
            2: dict(ref='TT',query='TT'),
            4: dict(ref='AA',query='A='),
            6: dict(ref='C',query='C'),
            7: dict(ref='A',query='A'),
            8: dict(ref='T',query='T'),
            9: dict(ref='A',query='A'),
            10: dict(ref='G',query='-')
            }

    def test07(self):

        #             01234567890
        self.ref =   'GCTTCAGA=C'
        self.query = '-ATT---AA-'

        self.expected = {
            0: dict(ref='G',query='-'),
            1: dict(ref='C',query='A'),
            2: dict(ref='TT',query='TT'),
            4: dict(ref='C',query='-'),
            5: dict(ref='A',query='-'),
            6: dict(ref='G',query='-'),
            7: dict(ref='A=',query='AA'),
            8: dict(ref='C', query='-')
            }


