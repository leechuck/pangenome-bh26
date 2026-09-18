import unittest
from evaluate_hgsvc import score as original_score
from evaluate_hgsvc_corrected import score,truth_slots,METHODS
from hgsvc_truth import GENES
from test_hgsvc_evaluation import calls


class CorrectedEvaluationTest(unittest.TestCase):
    def test_perfect_and_wrong_two_field_calls(self):
        cohort=[dict(donor='D',family='F',stratum='EAS')]
        truth={('D',g):[dict(status='exact_genomic',exact_genomic_alleles=[g+'*01:01:01:01'])]*2 for g in GENES}
        results={('D',m):dict(state='complete',calls=calls('01:01:01:01')) for m in METHODS}
        old_rows,_=original_score(cohort,results,truth,repetitions=100)
        self.assertEqual(sum(r['correct'] for r in old_rows if r['fields']==2),0)
        rows,_=score(cohort,results,truth,repetitions=100)
        self.assertTrue(all(r['correct']==1 for r in rows))
        self.assertEqual([r for r in rows if r['fields']==4],[r for r in old_rows if r['fields']==4])
        results['D','T1K']['calls']=calls('02:01:01:01')
        rows,_=score(cohort,results,truth,repetitions=100)
        self.assertTrue(all(r['correct']==0 for r in rows if r['method']=='T1K'))

    def test_eligibility_and_ambiguity_unchanged(self):
        records=[dict(status='exact_genomic',exact_genomic_alleles=['A*01:01:01:01','A*02:01:01:01'])]*2
        self.assertEqual(truth_slots(records,2),(frozenset(['01:01','02:01']),)*2)
        records[0]=dict(status='missing',exact_genomic_alleles=[])
        self.assertIsNone(truth_slots(records,2))


if __name__=='__main__':unittest.main()
