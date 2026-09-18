import unittest
from score_gourraud import family_interval

class FamilyIntervalTests(unittest.TestCase):
    def test_final_failure_is_bound_to_failed_attempt(self):
        import json,tempfile
        from pathlib import Path
        from score_series import terminal_failure,sha
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);m=p/'manifest.json';m.write_text(json.dumps({'status':'failed'}))
            self.assertFalse(terminal_failure(p))
            (p/'TERMINAL_FAILURE.json').write_text(json.dumps({'manifest_sha256':sha(m)}))
            self.assertTrue(terminal_failure(p))
            m.write_text(json.dumps({'status':'running'}))
            self.assertFalse(terminal_failure(p))
    def test_genotype_weighting_is_preserved(self):
        result=family_interval({'family_of_three':(3,15),'single':(-1,5)},reps=100)
        self.assertEqual(result['difference'],.1)
        self.assertEqual(result['families'],2)
        self.assertGreaterEqual(result['lower_95'],-.2)
        self.assertLessEqual(result['upper_95'],.2)
    def test_paired_zero_and_no_eligible_truth(self):
        result=family_interval({'a':(0,10),'b':(0,5)},reps=100)
        self.assertEqual((result['difference'],result['lower_95'],result['upper_95']),(0,0,0))
        self.assertIsNone(family_interval({}))
        self.assertIsNone(family_interval({'a':(0,0)}))

if __name__=='__main__':unittest.main()
