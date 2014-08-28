"""
Test sequtils module.
"""

import cPickle
import csv
import logging
import operator
import pprint
import sys

from bz2 import BZ2File
from collections import Counter
from os import path

from bioy_pkg import sequtils

from __init__ import TestBase, datadir as datadir

log = logging.getLogger(__name__)

sequtilsdir = path.join(datadir, 'sequtils')


class TestRunSsearch(TestBase):

    def test01(self):
        """
        Provide an output name
        """

        q, t = self.data('two.fasta'), self.data('ten.fasta')
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

        q, t = self.data('two.fasta'), self.data('ten.fasta')
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

        q, t = self.data('two.fasta'), self.data('ten.fasta')
        with sequtils.run_ssearch(q, t) as aligns:
            parsed = sequtils.parse_ssearch36(aligns)
            self.assertEqual(set(['H59735', 'T70875']),
                             {d['q_name'] for d in parsed})

        # should not still exist
        self.assertFalse(path.exists(aligns.name))

    def test04(self):
        with BZ2File(self.data('rle_100_left.ssearch.bz2')) as f:
            aligns = list(sequtils.parse_ssearch36(f))

            # 100 total sequences
            self.assertEqual(len(set(a['q_name'] for a in aligns)), 100)

            # searched against 4 primers
            self.assertEqual(len(aligns), 400)


class TestAllPairwise(TestBase):

    def test01(self):
        with open(self.data('five.fasta')) as f:
            seqs = list(sequtils.fastalite(f))
            pairs = list(sequtils.all_pairwise(seqs))
            self.assertEqual(len(pairs), (len(seqs) * (len(seqs) - 1)) / 2)
            self.assertEqual(
                [s.id for s in seqs], list(sequtils.names_from_pairs(pairs)))

    def test02(self):
        with open(self.data('two.fasta')) as f:
            seqs = list(sequtils.fastalite(f))
            pairs = list(sequtils.all_pairwise(seqs))
            self.assertEqual(len(pairs), (len(seqs) * (len(seqs) - 1)) / 2)
            self.assertEqual(
                [s.id for s in seqs], list(sequtils.names_from_pairs(pairs)))


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
        with sequtils.fasta_tempfile(self.seqs, dir=self.outdir) as f:
            self.assertTrue(path.exists(f))
        self.assertFalse(path.exists(f))


class TestEncodeAndDecode(TestBase):

    def test01(self):
        seq = 'TCTGGACCGTGTCTTTCAGTTCCAAAGTGTGACTGATCCATCCTCTCAGACC'

        e, c = sequtils.homoencode(seq)

        self.assertEquals(seq, sequtils.homodecode(e, c))


class TestRle(TestBase):

    def test01(self):
        rle, counts = sequtils.homoencode('AAA')
        self.assertEquals(rle, 'A')
        self.assertEquals(counts, [3])

    def test02(self):
        rle, counts = sequtils.homoencode('GCTTCAAACATA')
        self.assertEquals(rle, 'GCTCACATA')
        self.assertEquals(counts, [1, 1, 2, 1, 3, 1, 1, 1, 1])


class TestDecodeAlignment(TestBase):

    def test01(self):
        pass


class TestErrorCounting(TestBase):

    def tearDown(self):

        # print '\nref:  ', self.ref
        # print 'query:', self.query
        errors = list(sequtils.itemize_errors(self.ref, self.query))
        log.debug(sequtils.show_errors(errors))

        log.debug('\n' + pprint.pformat(self.expected))
        log.debug('\n' + pprint.pformat(errors))

        # indices into ref are as expected
        self.assertEquals(set(e['i']
                              for e in errors), set(self.expected.keys()))

        # individual positions are as expected
        for pos in errors:
            self.assertEquals(pos['ref'], self.expected[pos['i']]['ref'])
            self.assertEquals(pos['query'], self.expected[pos['i']]['query'])

    def test02(self):

        #             01234-56789
        self.ref = 'GATTA-CATA-'
        self.query = 'GCTTAACATAG'

        self.expected = {
            0: dict(ref='G', query='G'),
            1: dict(ref='A', query='C'),
            2: dict(ref='TT', query='TT'),
            4: dict(ref='A-', query='AA'),
            5: dict(ref='C', query='C'),
            6: dict(ref='A', query='A'),
            7: dict(ref='T', query='T'),
            8: dict(ref='A', query='A')
        }

    def test03(self):

        #             01234-56789
        self.ref = 'GATT-ACATA-'
        self.query = 'GCTTAACATAG'

        self.expected = {
            0: dict(ref='G', query='G'),
            1: dict(ref='A', query='C'),
            2: dict(ref='TT', query='TT'),
            4: dict(ref='-A', query='AA'),
            5: dict(ref='C', query='C'),
            6: dict(ref='A', query='A'),
            7: dict(ref='T', query='T'),
            8: dict(ref='A', query='A')
        }

    def test04(self):

        #             01234567890
        self.ref = 'GCTTAACATAG'
        self.query = 'GATT-ACATA-'

        self.expected = {
            0: dict(ref='G', query='G'),
            1: dict(ref='C', query='A'),
            2: dict(ref='TT', query='TT'),
            4: dict(ref='AA', query='-A'),
            6: dict(ref='C', query='C'),
            7: dict(ref='A', query='A'),
            8: dict(ref='T', query='T'),
            9: dict(ref='A', query='A'),
            10: dict(ref='G', query='-')
        }

    def test05(self):

        #              0123456789
        self.ref = '-CTTAACATAG'
        self.query = 'GATT-ACATAG'

        self.expected = {
            0: dict(ref='C', query='A'),
            1: dict(ref='TT', query='TT'),
            3: dict(ref='AA', query='-A'),
            5: dict(ref='C', query='C'),
            6: dict(ref='A', query='A'),
            7: dict(ref='T', query='T'),
            8: dict(ref='A', query='A'),
            9: dict(ref='G', query='G')
        }

    def test06(self):

        #             01234567890
        self.ref = 'GCTTAACATAG'
        self.query = '-ATTA=CATA-'

        self.expected = {
            0: dict(ref='G', query='-'),
            1: dict(ref='C', query='A'),
            2: dict(ref='TT', query='TT'),
            4: dict(ref='AA', query='A='),
            6: dict(ref='C', query='C'),
            7: dict(ref='A', query='A'),
            8: dict(ref='T', query='T'),
            9: dict(ref='A', query='A'),
            10: dict(ref='G', query='-')
        }

    def test07(self):

        #             01234567890
        self.ref = 'GCTTCAGA=C'
        self.query = '-ATT---AA-'

        self.expected = {
            0: dict(ref='G', query='-'),
            1: dict(ref='C', query='A'),
            2: dict(ref='TT', query='TT'),
            4: dict(ref='C', query='-'),
            5: dict(ref='A', query='-'),
            6: dict(ref='G', query='-'),
            7: dict(ref='A=', query='AA'),
            8: dict(ref='C', query='-')
        }


class TestAsciiEncoding(TestBase):

    def test01(self):
        v = range(79)
        self.assertEqual(v, sequtils.from_ascii(sequtils.to_ascii(v)))

    def test02(self):
        v = range(80)
        self.assertRaises(ValueError, sequtils.to_ascii, v)


class TestFastaLite(TestBase):

    def test01(self):
        with open(
            self.data('five.fasta')) as f, open(
                self.data('five.fasta')) as r:
            seqs = sequtils.fastalite(f)
            raw = r.read()
            fasta = ''
            for seq in seqs:
                fasta += '>{}\n{}\n'.format(seq.description, seq.seq)
                log.debug('{}'.format(seq))

        self.assertEquals(
            ''.join(raw).replace('\n', ''), fasta.replace('\n', ''))

    def test02(self):
        with open(self.data('five.fasta')) as f:
            seqs = sequtils.fastalite(f)
            for seq in seqs:
                pass


class TestParseClusters(TestBase):

    def test01(self):
        infile = self.data('clusters.uc')
        with open(infile) as f:
            cluster_ids, cluster_sizes = sequtils.parse_uc(f)

        counter = Counter()
        for cluster, count in cluster_sizes.items():
            counter[count] += 1

        # most of the clusters are singletons
        self.assertEquals(counter.most_common(1)[0][0], 1)


class TestCompoundAssignment(TestBase):

    thisdatadir = path.join(sequtilsdir, 'TestCompoundAssignment')

    assignments = path.join(thisdatadir, 'assignments.pkl.bz2')
    assignments = BZ2File(assignments)
    assignments = cPickle.load(assignments)

    taxonomy = BZ2File(path.join(datadir, 'taxonomy.csv.bz2'))
    taxonomy = {t['tax_id']: t for t in csv.DictReader(taxonomy)}

    def test01(self):
        """
        test basic behavior
        """

        taxonomy = self.taxonomy
        thisdatadir = self.thisdatadir

        this_test = sys._getframe().f_code.co_name

        compound_assignments_ref = path.join(thisdatadir,
                                             this_test,
                                             'compound_names.pkl.bz2')
        compound_assignments_ref = BZ2File(compound_assignments_ref)
        compound_assignments_ref = cPickle.load(compound_assignments_ref)

        compound = lambda x: sequtils.compound_assignment(x, taxonomy)
        compound_assignments = map(compound, self.assignments)

        self.assertEquals(compound_assignments, compound_assignments_ref)

    def test02(self):
        """
        test no data
        """

        taxonomy = self.taxonomy

        compound = lambda x: sequtils.compound_assignment(x, taxonomy)
        compound_assignments = map(compound, [])

        self.assertEquals(compound_assignments, [])

    def test03(self):
        """
        test no tax info
        """

        for a in self.assignments:
            self.assertRaises(TypeError, sequtils.compound_assignment, a, {})


class TestCondenseAssignment(TestBase):

    thisdatadir = path.join(sequtilsdir, 'TestCondenseAssignment')

    assignments = path.join(thisdatadir, 'assignments.pkl.bz2')
    assignments = BZ2File(assignments)
    assignments = cPickle.load(assignments)

    taxonomy = BZ2File(path.join(datadir, 'taxonomy.csv.bz2'))
    taxonomy = {t['tax_id']: t for t in csv.DictReader(taxonomy)}

    def test01(self):
        """
        test basic behavior, max_size = 3 behavior
        """

        taxonomy = self.taxonomy
        thisdatadir = self.thisdatadir

        this_test = sys._getframe().f_code.co_name

        condensed_assignments_ref = path.join(thisdatadir,
                                              this_test,
                                              'assignments.pkl.bz2')
        condensed_assignments_ref = BZ2File(condensed_assignments_ref)
        condensed_assignments_ref = cPickle.load(condensed_assignments_ref)

        condense_assignments = lambda x: sequtils.condense_ids(
            x, taxonomy, max_size=3)
        condensed_assignments = map(condense_assignments, self.assignments)

        self.assertEquals(condensed_assignments, condensed_assignments_ref)

    def test02(self):
        """
        test max_size = 1 behavior
        """

        taxonomy = self.taxonomy
        thisdatadir = self.thisdatadir

        this_test = sys._getframe().f_code.co_name

        condensed_assignments_ref = path.join(thisdatadir,
                                              this_test,
                                              'assignments.pkl.bz2')
        condensed_assignments_ref = BZ2File(condensed_assignments_ref)
        condensed_assignments_ref = cPickle.load(condensed_assignments_ref)

        condense_assignments = lambda x: sequtils.condense_ids(
            x, taxonomy, max_size=1)
        condensed_assignments = map(condense_assignments, self.assignments)

        self.assertEquals(condensed_assignments, condensed_assignments_ref)

    def test03(self):
        """
        test max_size = 0
        """

        taxonomy = self.taxonomy
        thisdatadir = self.thisdatadir

        this_test = sys._getframe().f_code.co_name

        condensed_assignments_ref = path.join(thisdatadir,
                                              this_test,
                                              'assignments.pkl.bz2')
        condensed_assignments_ref = BZ2File(condensed_assignments_ref)
        condensed_assignments_ref = cPickle.load(condensed_assignments_ref)

        condense_assignments = lambda x: sequtils.condense_ids(
            x, taxonomy, max_size=0)
        condensed_assignments = map(condense_assignments, self.assignments)

        self.assertEquals(condensed_assignments, condensed_assignments_ref)

    def test04(self):
        """
        test no taxonomy
        """

        for a in self.assignments:
            self.assertRaises(TypeError, sequtils.condense_ids, a, None)

    def test05(self):
        """
        test wrong floor_rank
        """

        taxonomy = self.taxonomy

        for a in self.assignments:
            self.assertRaises(
                TypeError,
                sequtils.condense_ids,
                a,
                taxonomy,
                floor_rank='fake_rank')

    def test06(self):
        """
        test wrong ceiling
        """

        taxonomy = self.taxonomy

        for a in self.assignments:
            self.assertRaises(
                TypeError,
                sequtils.condense_ids,
                a,
                taxonomy,
                ceiling_rank='fake_rank')

    def test07(self):
        """
        test ceiling too high
        """

        taxonomy = self.taxonomy

        for a in self.assignments:
            self.assertRaises(TypeError,
                              sequtils.condense_ids,
                              a,
                              taxonomy,
                              ceiling_rank='species',
                              floor_rank='superkingdom')


class TestCorrectCopyNumbers(TestBase):

    thisdatadir = path.join(sequtilsdir, 'TestCorrectCopyNumbers')

    assignments = path.join(thisdatadir, 'assignments.pkl.bz2')
    assignments = BZ2File(assignments)
    assignments = cPickle.load(assignments)

    copy_numbers = BZ2File(path.join(datadir, 'rrnDB_16S_copy_num.csv.bz2'))
    copy_numbers = {c['tax_id']: c['median']
                    for c in csv.DictReader(copy_numbers)}

    def test01(self):
        """
        test basic behavior
        """

        thisdatadir = self.thisdatadir

        this_test = sys._getframe().f_code.co_name

        corrections_ref = path.join(thisdatadir,
                                    this_test,
                                    'copy_numbers.pkl.bz2')
        corrections_ref = BZ2File(corrections_ref)
        corrections_ref = cPickle.load(corrections_ref)

        corrections = lambda x: sequtils.correct_copy_numbers(
            x, self.copy_numbers)
        corrections = map(corrections, self.assignments)

        self.assertEquals(corrections_ref, corrections)

    def test02(self):
        """
        test empty copy_corrections
        """

        thisdatadir = self.thisdatadir

        this_test = sys._getframe().f_code.co_name

        corrections_ref = path.join(thisdatadir,
                                    this_test,
                                    'copy_numbers.pkl.bz2')
        corrections_ref = BZ2File(corrections_ref)
        corrections_ref = cPickle.load(corrections_ref)

        corrections = lambda x: sequtils.correct_copy_numbers(
            x, self.copy_numbers)
        corrections = map(corrections, self.assignments)

        self.assertEquals(corrections_ref, corrections)

    def test03(self):
        """
        test empty assignments
        """

        self.assertRaises(
            TypeError, sequtils.correct_copy_numbers, [], self.copy_numbers)
