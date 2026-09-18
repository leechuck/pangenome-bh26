import unittest
import numpy as np
from joint_model import fit, emission_matrix


class JointModelTests(unittest.TestCase):
    def test_heterozygous_pair(self):
        result=fit([[1,0],[1,0],[0,1],[0,1]],['A','A'])
        self.assertTrue(result['converged'])
        self.assertEqual(result['genes']['A']['pairs'],[[0,1]])

    def test_homozygous_pair(self):
        self.assertEqual(fit([[1,.01]]*10,['A','A'])['genes']['A']['pairs'],[[0,0]])

    def test_identical_candidates_remain_ambiguous(self):
        self.assertEqual(len(fit([[1,1]]*5,['A','A'])['genes']['A']['pairs']),3)

    def test_competing_genes_share_fragment_mass(self):
        result=fit([[1,0]]*6+[[0,1]]*2+[[1,1]]*4,['A','B'])
        self.assertAlmostEqual(result['genes']['A']['gene_fraction'],.75,places=4)
        self.assertAlmostEqual(sum(x['gene_fraction'] for x in result['genes'].values()),1)
        self.assertTrue(all(b>=a-1e-8 for a,b in zip(result['objective_trace'],result['objective_trace'][1:])))

    def test_unassigned_fragments_cannot_create_gene_call(self):
        result=fit([[0,0],[0,0]],['A','B'])
        self.assertEqual(result['genes'],{})
        self.assertEqual(result['unassigned_fragments'],2)

    def test_duplicate_placements_do_not_change_emission(self):
        fragment={'candidates':{'p':[{'alignment_score':160}],'q':[{'alignment_score':150}]}}
        before=emission_matrix([fragment],['p','q'],{'p':100,'q':100})
        fragment['candidates']['p']*=10
        np.testing.assert_array_equal(before,emission_matrix([fragment],['p','q'],{'p':100,'q':100}))


if __name__=='__main__':
    unittest.main()
