import unittest
from run_graph_development_batch import task_coordinates,CONDITIONS


class BatchCoordinatesTests(unittest.TestCase):
    def test_bootstrap_is_additive_unsampled_for_every_gene(self):
        for index in range(56):
            donor,condition,gene=task_coordinates('bootstrap',index,7)
            self.assertEqual((donor,gene),(index//8,index%8))
            self.assertEqual(CONDITIONS[condition],('hprc_asian','unsampled'))

    def test_all_four_conditions_all_genes_exactly_once(self):
        actual={task_coordinates('mapping',index,7) for index in range(224)}
        self.assertEqual(actual,{(d,c,g) for d in range(7) for c in range(4) for g in range(8)})
        for index in range(28):
            self.assertEqual(task_coordinates('join',index,7)[:2],divmod(index,4))

    def test_out_of_bounds_rejected(self):
        for stage,maximum in [('coverage',7),('bootstrap',56),('mapping',224),('join',28)]:
            for index in (-1,maximum):
                with self.assertRaises(ValueError):task_coordinates(stage,index,7)
