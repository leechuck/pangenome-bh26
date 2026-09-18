import unittest
from graph_pair_refinement_native_locus import assign_fragment,native_locus_allowed


class NativeLocusTests(unittest.TestCase):
    def test_exclusive_paralog_assignment_vetoes_apparent_graph_evidence(self):
        record={'native_candidates':{'HLA-E*01:01':[{'adjusted_weight':1.}]},
                'graph_sources':{'A':{'candidates':{'p':[{'alignment_score':300}]}}}}
        paths={'p':dict(gene='A',effective_length=7000,candidates=['A*68:01:02:01'])}
        self.assertIsNone(assign_fragment(record,paths))
        record['native_candidates']={'HLA-A*68:01':[{'adjusted_weight':1.}]}
        self.assertEqual(assign_fragment(record,paths)[0],'A')

    def test_no_positive_native_evidence_does_not_discard_graph_only_reads(self):
        for native in ({},{'HLA-E*01:01':[{'adjusted_weight':0.}]}):
            self.assertTrue(native_locus_allowed({'native_candidates':native},'A'))

    def test_mixed_native_loci_are_not_misrepresented_as_exclusive(self):
        native={a:[{'adjusted_weight':.5}] for a in ('HLA-A*68:01','HLA-E*01:01')}
        self.assertTrue(native_locus_allowed({'native_candidates':native},'A'))
        self.assertFalse(native_locus_allowed({'native_candidates':native},'C'))
