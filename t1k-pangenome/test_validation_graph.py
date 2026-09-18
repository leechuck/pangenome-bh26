import unittest
from run_validation_graph import coordinates,STAGES,GENES


class ValidationGraphTests(unittest.TestCase):
    def test_all_donors_loci_and_panels_once(self):
        for stage in STAGES:
            width=16 if stage in ('mapping','evidence') else 8 if stage in ('bootstrap','bootstrap_evidence') else 2 if stage in ('join','refine') else 1
            rows=[coordinates(stage,i,84) for i in range(84*width)]
            self.assertEqual(len(set(rows)),84*width)
            self.assertEqual({r[0] for r in rows},set(range(84)))
            if width==16:
                self.assertEqual(set(rows),{(d,p,g) for d in range(84) for p in ('hprc','hprc_asian') for g in GENES})
            if width==8:self.assertEqual({r[1] for r in rows},{'hprc_asian'})
            for index in (-1,84*width):
                with self.assertRaises(ValueError):coordinates(stage,index,84)
