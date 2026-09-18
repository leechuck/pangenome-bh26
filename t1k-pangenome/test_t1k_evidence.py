import unittest
from t1k_evidence import parse


class EvidenceTests(unittest.TestCase):
    def test_competing_genes_and_native_weights_retained(self):
        rows = list(parse(['r\tHLA-A*01:01\t10\t100\t1\t1\t0.25\n',
                           'r\tHLA-B*01:01\t20\t110\t0.5\t0.001\t0.125\n']))
        self.assertEqual(len(rows),1)
        self.assertEqual(len(rows[0]['candidates']),2)
        self.assertEqual(rows[0]['candidates']['HLA-B*01:01'][0]['quality'],.001)

    def test_identical_assignment_not_counted_twice(self):
        line = 'r\tHLA-A*01:01\t10\t100\t1\t1\t1\n'
        self.assertEqual(len(list(parse([line,line]))[0]['candidates']['HLA-A*01:01']),1)

    def test_noncontiguous_fragment_is_rejected(self):
        line = '{}\tHLA-A*01:01\t10\t100\t1\t1\t1\n'
        with self.assertRaises(ValueError):
            list(parse([line.format(x) for x in ('r','s','r')]))

    def test_standard_export_and_nonfinite_weights_rejected(self):
        for line in ('r\tHLA-A*01:01\t10\t100\n',
                     'r\tHLA-A*01:01\t10\t100\tnan\t1\t1\n'):
            with self.assertRaises(ValueError):
                list(parse([line]))


if __name__ == '__main__':
    unittest.main()
