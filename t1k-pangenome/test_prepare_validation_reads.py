import io
import unittest
from prepare_validation_reads import verify_pairs


class PairVerificationTests(unittest.TestCase):
    def test_mate_suffixes(self):
        self.assertEqual(verify_pairs(io.StringIO('@r/1\nAC\n+\nII\n'),
                                      io.StringIO('@r/2\nGT\n+\nII\n')),1)

    def test_truncated_quality(self):
        with self.assertRaises(ValueError):
            verify_pairs(io.StringIO('@r\nAC\n+\nI\n'),io.StringIO('@r\nGT\n+\nII\n'))

    def test_missing_mate(self):
        with self.assertRaises(ValueError):
            verify_pairs(io.StringIO('@r\nAC\n+\nII\n'),io.StringIO(''))

    def test_reordered_mates(self):
        with self.assertRaises(ValueError):
            verify_pairs(io.StringIO('@r\nAC\n+\nII\n'),io.StringIO('@s\nGT\n+\nII\n'))


if __name__ == '__main__':
    unittest.main()
