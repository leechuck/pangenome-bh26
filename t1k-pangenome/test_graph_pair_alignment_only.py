import unittest
from graph_pair_refinement import assign_fragment as frozen_assign
from graph_pair_refinement_alignment_only import assign_fragment


class AlignmentOnlyTests(unittest.TestCase):
    def test_identical_alignment_does_not_gain_support_from_path_length(self):
        paths={p:dict(gene='C',effective_length=length,candidates=[p])
               for p,length in [('long',8000),('short',4000)]}
        record={'graph_sources':{'C':{'candidates':{
            p:[{'alignment_score':300}] for p in paths}}}}
        self.assertEqual(frozen_assign(record,paths)[1],{'long':.5,'short':1.})
        self.assertEqual(assign_fragment(record,paths)[1],{'long':1.,'short':1.})
        self.assertEqual(paths['long']['effective_length'],8000)

    def test_locus_competition_still_rejects_ambiguous_fragment(self):
        paths={p:dict(gene=p,effective_length=length,candidates=[p])
               for p,length in [('A',8000),('C',4000)]}
        record={'graph_sources':{'both':{'candidates':{
            'A':[{'alignment_score':300}],'C':[{'alignment_score':295}]}}}}
        self.assertIsNone(assign_fragment(record,paths))
