import unittest
from audit_reads import resolve


class ResolveReadsTest(unittest.TestCase):
    def test_local_path_resolves_through_public_manifest(self):
        self.assertEqual(resolve('HG1', '/ddbj/HG1.cram', {'HG1':'data/HG1.final.cram'}),
                         'https://s3.amazonaws.com/1000genomes/data/HG1.final.cram')

    def test_conflicting_recorded_source_rejected(self):
        with self.assertRaises(ValueError):
            resolve('HG1', 'https://different/HG1.cram', {'HG1':'data/HG1.final.cram'})

    def test_wrong_donor_path_rejected(self):
        with self.assertRaises(ValueError):
            resolve('HG1', '/ddbj/HG1.cram', {'HG1':'data/HG2.final.cram'})


if __name__ == '__main__':
    unittest.main()
