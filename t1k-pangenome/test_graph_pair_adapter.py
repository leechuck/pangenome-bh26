import hashlib
import unittest
from run_graph_pair_refinement import candidate_metadata,native_pair


def record(path,genomic,cds):
    return dict(path=path,sequence_sha256=hashlib.sha256(('A'*500).encode()).hexdigest(),genomic_labels=genomic,cds_labels=cds)


class AdapterTests(unittest.TestCase):
    def test_context_multiplicity_does_not_multiply_alleles(self):
        allele='A*01:01:01:01'
        records=[record(p,[allele],[allele]) for p in ('p','q')]
        paths,candidates=candidate_metadata('A',records,{'p':'A'*500,'q':'A'*500},[allele,allele],200)
        self.assertEqual(len(candidates),1)
        self.assertEqual(paths['p']['candidates'],paths['q']['candidates'])

    def test_cds_only_and_incompatible_paths_remain_nuisance(self):
        records=[record('p',[],['A*01:01:01:01']),record('q',['A*02:01:01:01'],['A*02:01:01:01'])]
        paths,candidates=candidate_metadata('A',records,{'p':'A'*500,'q':'A'*500},['A*01:01:01:01']*2,200)
        self.assertEqual(len(paths),2)
        self.assertTrue(all(c['label'] is None and c['families']==['A*01:01'] for c in candidates))

    def test_changed_sequence_rejected(self):
        with self.assertRaises(ValueError):
            candidate_metadata('A',[record('p',['A*01:01:01:01'],[])],{'p':'C'*500},['A*01:01:01:01']*2,200)

    def test_ambiguous_native_slot_remains_unresolved(self):
        call=dict(resolutions={'4':[dict(unresolved=False,alleles=['A*01:01:01:01','A*01:01:01:02'])]*2})
        self.assertEqual(native_pair(call),[])
