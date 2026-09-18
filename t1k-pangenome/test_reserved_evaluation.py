import unittest
from evaluate_reserved import METHODS,score


class Names:
    def g_groups(self,*args):return {'synthetic'}


class ReservedEvaluationTests(unittest.TestCase):
    def setUp(self):
        self.cohort=[dict(donor=d,family=d,stratum='EAS') for d in ('one','two')]
        allele='A*01:01:01:01'
        self.truth={(d,'A'):[dict(exact_genomic_alleles=allele,exact_cds_alleles=allele)]*2 for d in ('one','two')}
        call=dict(resolutions={str(n):[dict(alleles=[':'.join(allele.split(':')[:n])],unresolved=False)]*2 for n in (2,4)})
        self.results={(d,m):dict(state='complete',calls={'A':call}) for d in ('one','two') for m in METHODS}

    def test_failed_donor_stays_in_denominator_and_secondary_contrasts_separate(self):
        self.results['two','graph_hprc_asian']=dict(state='failed',calls={})
        rows,report=score(self.cohort,self.results,self.truth,Names(),100)
        r=next(r for r in report['summary'] if r['method']=='graph_hprc_asian' and r['fields']==4)
        self.assertEqual((r['eligible'],r['correct'],r['failed_donors']),(2,1,1))
        self.assertEqual(len(rows),2*4*8*2)
        self.assertNotIn('protocol_pass',report['contrasts']['graph_hprc'])
        self.assertFalse(report['contrasts']['T1K']['protocol_pass'])

    def test_pending_or_missing_results_cannot_be_silently_omitted(self):
        self.results['two','graph_hprc_asian']['state']='running'
        with self.assertRaises(ValueError):score(self.cohort,self.results,self.truth,Names(),100)
        self.results.pop(('two','graph_hprc_asian'))
        with self.assertRaises(ValueError):score(self.cohort,self.results,self.truth,Names(),100)
