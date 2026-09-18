from pathlib import Path
import tempfile
import unittest
from projection_graph import normalize
from path_support import PathIndex


class ProjectionGraphTests(unittest.TestCase):
    def test_split_nodes_and_reverse_path_preserved(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);source=root/'export.gfa';reference=root/'reference.fa';output=root/'graph.gfa'
            source.write_text('S\t1\tAAA\nS\t7\tCC\nP\tp#0\t1+,7-\t*\n')
            reference.write_text('>p\nAAAGG\n')
            report=normalize(source,output,reference)
            self.assertEqual(report['path_aliases'],{'p#0':'p'})
            index=PathIndex.from_gfa(output)
            placements,status=index.placements({'path':{'mapping':[
                {'position':{'node_id':7,'is_reverse':True},'edit':[{'from_length':2}]}]}})
            self.assertEqual(status,'supported')
            self.assertEqual(placements,[dict(path='p',reverse=False,start=3,end=5)])

    def test_sequence_change_rejected_before_output(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);source=root/'export.gfa';reference=root/'reference.fa';output=root/'graph.gfa'
            source.write_text('S\t1\tAAA\nP\tp#0\t1+\t*\n');reference.write_text('>p\nAAT\n')
            with self.assertRaises(ValueError):normalize(source,output,reference)
            self.assertFalse(output.exists())
