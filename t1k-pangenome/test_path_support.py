import unittest
from path_support import PathIndex


def alignment(nodes, lengths, offsets=None):
    return {'path':{'mapping':[{'position':{'node_id':abs(n),'is_reverse':n<0,
                  'offset':(offsets or [0]*len(nodes))[i]},'edit':[{'from_length':lengths[i]}]}
                  for i,n in enumerate(nodes)]}}


class PathSupportTests(unittest.TestCase):
    def setUp(self):
        self.index=PathIndex({1:10,2:5,3:7,4:5},{'a':[1,2,3],'b':[1,4,3]})

    def test_shared_sequence_retains_both_paths(self):
        placements,status=self.index.placements(alignment([1],[5],[2]))
        self.assertEqual({p['path'] for p in placements},{'a','b'})
        self.assertEqual(status,'supported')

    def test_variant_walk_selects_compatible_path(self):
        placements,_=self.index.placements(alignment([1,2],[3,4],[7,0]))
        self.assertEqual(placements,[dict(path='a',reverse=False,start=7,end=14)])

    def test_reverse_coordinates(self):
        placements,_=self.index.placements(alignment([-3,-2],[5,4],[2,0]))
        self.assertEqual(placements,[dict(path='a',reverse=True,start=11,end=20)])

    def test_recombination_is_not_forced_onto_path(self):
        self.assertEqual(self.index.placements(alignment([2,4],[5,5]))[1],'no_observed_path')

    def test_disconnected_partial_nodes_rejected(self):
        self.assertEqual(self.index.placements(alignment([1,2],[3,4]))[1],
                         'noncontiguous_mapping_bounds')

    def test_wrong_graph_rejected(self):
        with self.assertRaises(ValueError):
            self.index.placements(alignment([99],[5]))


if __name__ == '__main__':
    unittest.main()
