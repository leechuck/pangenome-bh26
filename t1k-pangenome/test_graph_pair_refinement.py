import unittest
from graph_pair_refinement import refine_gene,assign_fragment,internal_fragment


class GraphRefinementTests(unittest.TestCase):
    def test_internal_pair_rejects_flanking_mate_even_if_other_mate_overlaps(self):
        paths={'p':dict(gene='A',gene_span=[100,1000])}
        placement=dict(first_start=150,first_end=250,second_start=950,second_end=1050)
        record=dict(graph_sources={'g':dict(candidates={'p':[placement]})})
        self.assertFalse(internal_fragment(record,paths,'A'))
        placement.update(second_start=700,second_end=800)
        self.assertTrue(internal_fragment(record,paths,'A'))

    def test_boundary_ambiguity_rejects_whole_fragment_without_deleting_candidates(self):
        paths={p:dict(gene='A',gene_span=[100,1000]) for p in ('p','q')}
        inside=dict(first_start=150,first_end=250,second_start=700,second_end=800)
        outside=dict(inside,first_start=50)
        record=dict(graph_sources={'g':dict(candidates={'p':[inside],'q':[outside]})})
        self.assertFalse(internal_fragment(record,paths,'A'))
        self.assertEqual(set(record['graph_sources']['g']['candidates']),{'p','q'})
        paths['q']['gene']='B'
        self.assertTrue(internal_fragment(record,paths,'A'))

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
