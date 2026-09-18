import unittest
from library_calibration import candidate,select,estimate


class LibraryCalibrationTests(unittest.TestCase):
    def test_path_duplicates_do_not_add_fragments(self):
        p=dict(alignment_score=300,insert_size=400)
        fragment=dict(read_names=['r/1','r/2'],candidates={'a':[p],'b':[p]})
        selected,_=select([candidate(fragment,'A')])
        self.assertEqual(len(selected),1)
        self.assertEqual(selected[0]['span'],400)

    def test_close_competing_locus_and_ambiguous_span_rejected(self):
        rows=[dict(fragment='r',gene='A',score=300,span=400),
              dict(fragment='r',gene='B',score=295,span=450),
              dict(fragment='s',gene='A',score=300,span=None)]
        selected,counts=select(rows)
        self.assertFalse(selected)
        self.assertEqual(counts,dict(competing_locus=1,ambiguous_insert=1))

    def test_unique_best_locus_contributes_once(self):
        rows=[dict(fragment='r',gene='A',score=300,span=400),
              dict(fragment='r',gene='B',score=280,span=450)]
        self.assertEqual(select(rows)[0],[rows[0]])
        with self.assertRaises(ValueError):select([rows[0],rows[0]])

    def test_robust_estimate_and_insufficient_evidence(self):
        spans=list(range(300,501))*2
        fit=estimate(spans+[10000])
        self.assertTrue(fit['usable']);self.assertEqual(fit['fragment_mean'],400)
        self.assertEqual(fit['trimmed_fragments'],1)
        self.assertFalse(estimate([400]*300)['usable'])
        self.assertFalse(estimate(spans[:100])['usable'])
