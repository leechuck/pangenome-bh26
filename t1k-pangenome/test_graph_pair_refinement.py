import unittest
from graph_pair_refinement import refine_gene,assign_fragment


class GraphRefinementTests(unittest.TestCase):
    def setUp(self):
        self.a='DQA1*01:04:01:01';self.b='DQA1*01:04:01:02'
        self.candidates=[dict(label=a,families=['DQA1*01:04']) for a in (self.a,self.b)]

    def test_read_evidence_can_correct_false_heterozygosity(self):
        r=refine_gene([self.a,self.b],self.candidates,[[.001,1]]*100)
        self.assertTrue(r['changed']);self.assertEqual(r['pair'],[self.b,self.b])

    def test_missing_native_candidate_keeps_ipd_fallback(self):
        r=refine_gene([self.a,self.b],self.candidates[:1],[[1]]*100)
        self.assertFalse(r['changed']);self.assertIn('fallback',r['reason'])

    def test_ambiguity_and_unlabelled_winner_do_not_force_calls(self):
        self.assertFalse(refine_gene([self.a,self.b],self.candidates,[[1,1]]*100)['changed'])
        candidates=self.candidates+[dict(label=None,families=['DQA1*01:04'])]
        r=refine_gene([self.a,self.b],candidates,[[.001,.001,1]]*100)
        self.assertFalse(r['changed']);self.assertIn('unlabelled',r['reason'])

    def test_different_two_field_pair_cannot_be_substituted(self):
        c=self.candidates+[dict(label='DQA1*02:01:01:01',families=['DQA1*02:01'])]
        r=refine_gene([self.a,self.b],c,[[1,1,1]]*100)
        self.assertFalse(r['changed'])

    def test_path_duplicates_do_not_add_evidence_and_competing_locus_rejected(self):
        paths={p:dict(gene='A',effective_length=100,candidates=['a']) for p in ('p','duplicate')}
        r=dict(graph_sources={'g':dict(candidates={'p':[dict(alignment_score=300)]})})
        before=assign_fragment(r,paths)
        r['graph_sources']['g']['candidates']['duplicate']=[dict(alignment_score=300)]
        self.assertEqual(assign_fragment(r,paths),before)
        paths['b']=dict(gene='B',effective_length=100,candidates=['b'])
        r['graph_sources']['g']['candidates']['b']=[dict(alignment_score=295)]
        self.assertIsNone(assign_fragment(r,paths))
