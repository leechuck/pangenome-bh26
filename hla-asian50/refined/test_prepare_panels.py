import unittest
from prepare_panels import clean_gt,remap_site
class MissingTests(unittest.TestCase):
 def test_preserve_partial(self):
  aa=['A','T','C'];self.assertEqual(clean_gt('1|.',aa),(1,None));self.assertEqual(clean_gt('.|2',aa),(None,2))
 def test_unknown_not_reference(self):
  aa=['A','T','N'];self.assertEqual(clean_gt('2|1',aa),(None,1));self.assertEqual(clean_gt('0/1',aa),(None,None));self.assertEqual(clean_gt('0|1',aa,True),(None,None))
 def test_union_retains_known_and_excludes_test_only(self):
  aa=['A','T','C','G'];base=[(0,1)];added=[(None,2)]
  self.assertEqual(remap_site(aa,base),([0,1],['0|1']))
  self.assertEqual(remap_site(aa,base+added),([0,1,2],['0|1','.|2']))
 def test_no_known_alt(self):self.assertIsNone(remap_site(['A','T'],[(0,None),(None,None)]))
if __name__=='__main__':unittest.main()
