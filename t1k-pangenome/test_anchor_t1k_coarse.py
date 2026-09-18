import copy
import unittest
from anchor_t1k_coarse import anchor


def call(two,four):
    return dict(copy_count=2,resolutions={str(n):[dict(alleles=[a],unresolved=False,raw=a) for a in v]
        for n,v in ((2,two),(4,four))})


class AnchorTests(unittest.TestCase):
    def test_accepts_refinement_and_preserves_original_coarse_slot_order(self):
        original={'A':call(['A*02:01','A*01:01'],['A*02:01:01:01','A*01:01:01:01'])}
        candidate={'A':call(['A*01:01','A*02:01'],['A*01:01:01:02','A*02:01:01:02'])}
        before=copy.deepcopy(original)
        out,_=anchor(original,candidate)
        self.assertEqual(original,before)
        self.assertEqual(out['A']['resolutions']['2'],original['A']['resolutions']['2'])
        self.assertEqual([s['alleles'][0] for s in out['A']['resolutions']['4']],['A*02:01:01:02','A*01:01:01:02'])

    def test_coarse_conflict_or_ambiguity_retains_entire_original_call(self):
        original={'A':call(['A*02:01','A*01:01'],['A*02:01:01:01','A*01:01:01:01'])}
        candidate={'A':call(['A*03:01','A*01:01'],['A*03:01:01:01','A*01:01:01:02'])}
        self.assertEqual(anchor(original,candidate)[0],original)
        original['A']['resolutions']['2'][0]['unresolved']=True
        self.assertEqual(anchor(original,candidate)[0],original)

    def test_missing_candidate_does_not_remove_original_gene(self):
        original={'A':call(['A*02:01','A*01:01'],['A*02:01:01:01','A*01:01:01:01'])}
        self.assertEqual(anchor(original,{})[0],original)
