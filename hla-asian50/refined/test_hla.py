import unittest
from prepare_hla import labels
from score_hla import compare,unique
class HlaTests(unittest.TestCase):
 def test_unphased_exact_pair(self):
  self.assertEqual(compare(('02:01','01:01'),[{'01:01'},{'02:01'}]),(1,1,2))
 def test_truth_ambiguity_not_prediction_choice(self):
  self.assertEqual(compare(('01:01','02:01'),[{'01:01','01:02'},{'02:01'}]),(1,1,2))
  self.assertIsNone(unique('A*01:01,A*01:02'))
  self.assertEqual(unique('A*01:01:01,A*01:01:02'),'01:01')
 def test_no_call_counts_as_failure(self):
  self.assertEqual(compare(('01:01',None),[{'01:01'},{'02:01'}]),(0,0,1))
 def test_label_requires_complete_cds(self):
  self.assertFalse(labels({'cds_complete':'0','exact_cds_alleles':'A*01:01:01'}))
  self.assertEqual(labels({'cds_complete':'1','exact_cds_alleles':'A*01:01:01;A*01:01:02N'}),{'A*01:01'})
if __name__=='__main__':unittest.main()
