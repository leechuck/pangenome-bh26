import unittest
from build_panel import arms

class ArmTests(unittest.TestCase):
    def setUp(self):
        self.donors={('GRCh38','REF'):('REF',['GRCh38#0'])}
        for i in range(10):self.donors[(f'H{i}','HPRC')]=('HPRC-Rest',[f'H{i}#1',f'H{i}#2'])
        for i in range(6):self.donors[(f'A{i}','APR')]=('APR',[f'A{i}#1',f'A{i}#2'])
    def test_sizes_matched_and_reference_kept(self):
        A=arms(self.donors)
        self.assertEqual(len(A['full']),33)
        for arm in ['hprc','asian_matched','random_matched']:
            self.assertEqual(len(A[arm]),21,arm)
            self.assertIn('GRCh38#0',A[arm])
    def test_asian_matched_contains_all_non_hprc(self):
        A=arms(self.donors)
        self.assertTrue({f'A{i}#{h}' for i in range(6) for h in '12'}<=A['asian_matched'])
        self.assertFalse({f'A{i}#1' for i in range(6)}&A['hprc'])
    def test_reassembly_of_hprc_donor_not_in_hprc_arm(self):
        self.donors[('H0','JaSaPaGe-Japanese')]=('JaSaPaGe-Japanese',['H0_J#1','H0_J#2'])
        A=arms(self.donors)
        self.assertNotIn('H0_J#1',A['hprc'])
        self.assertIn('H0_J#1',A['asian_matched'])
    def test_whole_donors_and_deterministic(self):
        A=arms(self.donors);B=arms(self.donors)
        self.assertEqual(A,B)
        for h in A['random_matched']:
            if h!='GRCh38#0':self.assertIn(h[:-1]+('2' if h[-1]=='1' else '1'),A['random_matched'])

if __name__=='__main__':unittest.main()
