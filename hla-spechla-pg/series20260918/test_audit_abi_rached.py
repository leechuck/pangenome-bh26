import unittest
from audit_abi_rached import alternatives


class PublishedAlleleTests(unittest.TestCase):
    def test_shorthand_and_null_suffix(self):
        self.assertEqual(alternatives('35:01/40N/42'), (['35:01', '35:40N', '35:42'], 'ambiguous'))

    def test_multiple_first_fields(self):
        self.assertEqual(alternatives('54:01 55:01/02 56:01'),
                         (['54:01', '55:01', '55:02', '56:01'], 'ambiguous'))

    def test_spreadsheet_corruption_is_not_an_allele(self):
        self.assertEqual(alternatives('0,0840277777777778'), ([], 'invalid_source_value'))

    def test_annotation_and_missing_are_not_guessed(self):
        self.assertEqual(alternatives('02:02*'), ([], 'source_asterisk_unresolved'))
        for value in ('', 'None'):
            self.assertEqual(alternatives(value), ([], 'missing'))


if __name__ == '__main__':
    unittest.main()
