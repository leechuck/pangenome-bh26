import unittest
from evaluate_hgsvc import score, METHODS
from hgsvc_truth import GENES


def calls(name):
    return {g:{'resolutions':{str(n):[{'unresolved':False,'alleles':[g+'*'+':'.join(name.split(':')[:n])]}]*2 for n in (2,4)}} for g in GENES}


class EvaluationTest(unittest.TestCase):
    def test_complete_grid_and_failed_predictions(self):
        cohort=[dict(donor='D1',family='F1',stratum='EAS'),dict(donor='D2',family='F2',stratum='AFR')]
        truth={(d,g):[dict(status='exact_genomic',exact_genomic_alleles=[g+'*01:01:01:01'])]*2 for d in ('D1','D2') for g in GENES}
        results={(d,m):dict(state='complete',calls=calls('01:01:01:01')) for d in ('D1','D2') for m in METHODS}
        results['D1','T1K']['calls']=calls('01:01:01:02')
        rows,report=score(cohort,results,truth,repetitions=100)
        self.assertEqual(len(rows),2*8*4*2)
        t1k=[s for s in report['summary'] if s['method']=='T1K' and s['fields']==4][0]
        self.assertEqual((t1k['correct'],t1k['eligible']),(8,16))
        results['D2','graph_hprc_asian']=dict(state='failed',calls={})
        rows,report=score(cohort,results,truth,repetitions=100)
        candidate=[s for s in report['summary'] if s['method']=='graph_hprc_asian' and s['fields']==4][0]
        self.assertEqual((candidate['correct'],candidate['eligible'],candidate['failed_donors']),(8,16,1))
        results.pop(('D2','graph_hprc'))
        with self.assertRaises(ValueError):score(cohort,results,truth,repetitions=100)


if __name__=='__main__':unittest.main()
