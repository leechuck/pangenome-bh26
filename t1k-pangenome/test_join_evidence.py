from contextlib import closing
import sqlite3
import unittest
from join_evidence import fragment_id, merge


class JoinTests(unittest.TestCase):
    def test_mate_suffixes_match_t1k_normalization(self):
        self.assertEqual(fragment_id(['read/1','read/2']),'read')
        self.assertEqual(fragment_id(['read','read']),'read')
        with self.assertRaises(ValueError):
            fragment_id(['read/1','other/2'])

    def test_same_fragment_across_methods_and_loci_is_one_row(self):
        native=[dict(fragment='r',candidates={'A*01:01':[],'B*01:01':[]})]
        graph=[dict(read_names=['r/1','r/2'],candidates={'pathA':[]})]
        with closing(sqlite3.connect(':memory:')) as connection:
            rows=list(merge(connection,native,[graph,graph]))
        self.assertEqual(len(rows),1)
        self.assertEqual(len(rows[0]['graph_sources']),2)
        self.assertEqual(len(rows[0]['native_candidates']),2)

    def test_missing_evidence_is_retained(self):
        with closing(sqlite3.connect(':memory:')) as connection:
            rows=list(merge(connection,[dict(fragment='native',candidates={'A':[]})],
                            [[dict(read_names=['graph','graph'],candidates={})]]))
        self.assertEqual({r['fragment'] for r in rows},{'native','graph'})

    def test_normalization_collision_is_rejected(self):
        graph=[dict(read_names=['r','r'],candidates={}),dict(read_names=['r/1','r/2'],candidates={})]
        with closing(sqlite3.connect(':memory:')) as connection, self.assertRaises(ValueError):
            list(merge(connection,[],[graph]))


if __name__=='__main__':
    unittest.main()
