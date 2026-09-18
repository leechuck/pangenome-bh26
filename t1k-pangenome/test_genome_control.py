import unittest
from build_genome_control import records, combine, invalid_annotations


class GenomeControlTests(unittest.TestCase):
    def test_missing_genomic_candidate_is_retained_without_padding(self):
        original = records('>HLA-A*01 1 0 2\nAAA\n>HLA-B*01 1 0 2\nCCC\n')
        genomic = records('>HLA-A*01 1 2 4\nGGAAATT\n>HLA-C*01 1 0 2\nTTT\n')
        result,replaced = combine(original,genomic)
        self.assertEqual(set(result),set(original))
        self.assertEqual(result['HLA-B*01'],original['HLA-B*01'])
        self.assertEqual(result['HLA-A*01'],genomic['HLA-A*01'])
        self.assertEqual(replaced,['HLA-A*01'])

    def test_invalid_annotation_and_duplicate_names_rejected(self):
        for source in ('>A 1 0 5\nAAA\n','>A 1 0 2\nAAA\n>A 1 0 2\nAAA\n'):
            with self.assertRaises(ValueError):
                records(source)

    def test_legacy_annotation_exception_is_retained_and_reported(self):
        original = records('>A 1 0 5\nAAA\n',validate_coordinates=False)
        self.assertEqual(set(invalid_annotations(original)),{'A'})
        result,replaced = combine(original,{})
        self.assertEqual(result,original)
        self.assertEqual(replaced,[])


if __name__ == '__main__':
    unittest.main()
