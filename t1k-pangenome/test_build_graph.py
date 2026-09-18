import tempfile
import unittest
from pathlib import Path
from build_graph import sequences, verify_paths


class PathPreservationTest(unittest.TestCase):
    def test_changed_base_is_rejected(self):
        with self.assertRaises(ValueError):
            verify_paths({'p':'ACGT'}, {'p':'ACGA'})

    def test_relabelled_path_is_rejected(self):
        with self.assertRaises(ValueError):
            verify_paths({'p':'ACGT'}, {'q':'ACGT'})

    def test_dropped_identical_path_is_rejected(self):
        with self.assertRaises(ValueError):
            verify_paths({'p':'ACGT','q':'ACGT'}, {'p':'ACGT'})

    def test_multiline_fasta_and_order(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)/'paths.fa'
            p.write_text('>b\nac\ngt\n>a\nTT\n')
            verify_paths({'a':'TT','b':'ACGT'}, sequences(p))

    def test_duplicate_names_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            p = Path(folder)/'paths.fa'
            p.write_text('>a\nAC\n>a\nTT\n')
            with self.assertRaises(ValueError):
                sequences(p)


if __name__ == '__main__':
    unittest.main()
