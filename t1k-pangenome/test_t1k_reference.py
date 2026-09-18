import unittest
from t1k_reference import header, resolve_alias


class ReferenceAdapterTests(unittest.TestCase):
    def test_t1k_inclusive_coordinates(self):
        self.assertEqual(header('A*PGabc', [[2,5],[7,9]], 10), '>A*PGabc 2 2 4 7 8\n')

    def test_overlapping_intervals_rejected(self):
        with self.assertRaises(ValueError):
            header('A*PGabc', [[2,5],[4,7]], 10)

    def test_cds_identity_does_not_imply_four_fields(self):
        aliases={'x':dict(genomic_labels=[],cds_labels=['A*01:01:01:01','A*01:01:01:02'])}
        self.assertEqual(resolve_alias('x',aliases,2)['alleles'],['A*01:01'])
        self.assertTrue(resolve_alias('x',aliases,4)['unresolved'])

    def test_genomic_ambiguity_preserved(self):
        aliases={'x':dict(genomic_labels=['A*01:01:01:01','A*01:01:01:02'])}
        result=resolve_alias('x',aliases,4)
        self.assertEqual(len(result['alleles']),2)
        self.assertFalse(result['unresolved'])

    def test_short_label_never_padded(self):
        self.assertTrue(resolve_alias('x',{'x':dict(genomic_labels=['A*01:01'])},4)['unresolved'])


if __name__ == '__main__':
    unittest.main()
