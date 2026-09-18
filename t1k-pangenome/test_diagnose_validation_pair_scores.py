import unittest
from diagnose_validation_pair_scores import pair_log
from graph_pair_refinement import refine_gene


class ScoreReplayTests(unittest.TestCase):
    def test_reproduces_frozen_pair_gain_including_zero_emissions(self):
        native='A*01:01:01:01';proposal='A*01:01:01:02'
        candidates=[dict(label=a,families=['A*01:01']) for a in (native,proposal)]
        rows=[[0.,1.]]*30+[[.2,.8]]*10
        result=refine_gene([native,native],candidates,rows,minimum_gap=0)
        self.assertTrue(result['changed'])
        replay=sum(pair_log(dict(zip((native,proposal),row)),[proposal,proposal])-
                   pair_log(dict(zip((native,proposal),row)),[native,native]) for row in rows)
        self.assertAlmostEqual(replay,result['score_gain_over_native'],places=10)
