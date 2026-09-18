import unittest
from decode_t1k import decode, decode_token


class DecodeTests(unittest.TestCase):
    def test_prefixed_context_alias(self):
        aliases={'HLA-A*PGx':dict(genomic_labels=['A*01:01:01:01'],cds_labels=[])}
        self.assertEqual(decode_token('HLA-A*PGx',aliases,4)['alleles'],['A*01:01:01:01'])

    def test_duplicate_gene_namespace_rejected(self):
        with self.assertRaises(ValueError):
            decode('HLA-A\t0\nA\t0',{})

    def test_unknown_alternative_is_not_discarded(self):
        result=decode_token('A*01:01:01:01,A*PGunknown',{},4)
        self.assertEqual(result['alleles'],['A*01:01:01:01'])
        self.assertTrue(result['unresolved'])

    def test_homozygous_context_decodes_both_copies(self):
        aliases={'A*PGx':dict(genomic_labels=['A*01:01:01:01'],cds_labels=['A*01:01:01:01'])}
        result=decode('A\t1\tA*PGx\t100\t60\t.\t0\t0',aliases)['A']['resolutions']['4']
        self.assertEqual(result[0]['alleles'],result[1]['alleles'])
        self.assertFalse(result[0]['unresolved'])

    def test_zero_confidence_is_unresolved(self):
        result=decode('A\t1\tA*01:01:01:01\t100\t0\t.\t0\t0',{})['A']['resolutions']['4']
        self.assertTrue(all(x['unresolved'] for x in result))

    def test_cross_gene_context_rejected(self):
        with self.assertRaises(ValueError):
            decode('A\t1\tB*01:01:01:01\t100\t60\t.\t0\t0',{})


if __name__ == '__main__':
    unittest.main()
