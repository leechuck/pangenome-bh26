import unittest
from hgsvc_truth import GENES, derive, reference_index, truth_slots


class TruthTest(unittest.TestCase):
    def setUp(self):
        self.refs = {g:{} for g in GENES}
        self.refs['A'] = {'ACGT':['A*01:01:01:01','A*01:01:01:02N']}

    def record(self, hap, core='ACGT'):
        return ('HG99999.'+hap+'_HLA-A_DELIBERATELY_WRONG_LABEL', 'G'*1000+core+'T'*1000)

    def test_exact_whole_core_not_header_or_substring(self):
        rows = derive([self.record('1'),self.record('2')], ['HG99999'], self.refs)
        self.assertEqual(len(rows),8)
        self.assertEqual(truth_slots(rows['HG99999','A'],4),
                         (frozenset(['A*01:01:01:01','A*01:01:01:02']),)*2)
        rows = derive([self.record('1','AACGT'),self.record('2')], ['HG99999'], self.refs)
        self.assertIsNone(truth_slots(rows['HG99999','A'],4))

    def test_missing_duplicate_ambiguous_fail_closed(self):
        for records, reason in (([self.record('2')],'missing'),
                                ([self.record('1'),self.record('1'),self.record('2')],'duplicate_haplotype_gene'),
                                ([self.record('1','ACNT'),self.record('2')],'ambiguous_core')):
            pair = derive(records,['HG99999'],self.refs)['HG99999','A']
            self.assertEqual(pair[0]['status'],reason)
            self.assertIsNone(truth_slots(pair,2))

    def test_short_names_never_padded(self):
        self.refs['A']['ACGT'].append('A*01:01')
        pair = derive([self.record('1'),self.record('2')],['HG99999'],self.refs)['HG99999','A']
        self.assertIsNone(truth_slots(pair,4))
        self.assertEqual(truth_slots(pair,2),(frozenset(['A*01:01']),)*2)

    def test_reference_retains_all_exact_names(self):
        result = reference_index([('HLA:1 A*01:01:01:01','ACGT'),
                                  ('HLA:2 A*01:01:01:02','ACGT'),
                                  ('HLA:3 A*02:01:01:01','ACNT')], 'A')
        self.assertEqual(result,{'ACGT':['A*01:01:01:01','A*01:01:01:02']})


if __name__ == '__main__':unittest.main()
