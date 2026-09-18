import unittest
from kmer_coverage import kmers, conserved_markers, estimate


class CoverageTests(unittest.TestCase):
    def test_reverse_complements_count_identically(self):
        sequence = 'AACCGTGATCTAG'
        reverse = sequence.translate(str.maketrans('ACGT','TGCA'))[::-1]
        self.assertEqual(sorted(kmers(sequence,5)),sorted(kmers(reverse,5)))

    def test_unknown_bases_break_windows(self):
        self.assertEqual(list(kmers('AAANCCC',3)),list(kmers('AAA',3))+list(kmers('CCC',3)))

    def test_repeated_kmers_are_not_single_copy_markers(self):
        self.assertEqual(conserved_markers(['AAAAAAAA'],3,1),set())

    def test_conservation_does_not_use_path_multiplicity_as_read_depth(self):
        self.assertEqual(conserved_markers(['AACCG'],3,1),conserved_markers(['AACCG']*10,3,1))

    def test_zero_markers_remain_in_coverage_estimate(self):
        result = estimate([0]*60+[100]*40)
        self.assertEqual(result['median_kmer_depth'],0)
        self.assertFalse(result['usable'])
        self.assertIsNone(result['personalization_coverage'])

    def test_sufficient_shared_depth_can_parameterize_sampling(self):
        self.assertEqual(estimate([65]*100)['personalization_coverage'],65)


if __name__ == '__main__':
    unittest.main()
