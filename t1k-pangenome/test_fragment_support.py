import unittest
from fragment_support import collect


def record(name,mate,slot,path,start,end,reverse):
    return dict(name=name, fragment_next={'name':mate} if slot==1 else None,
                fragment_prev={'name':mate} if slot==2 else None, score=160,
                placements=[dict(path=path,start=start,end=end,reverse=reverse)])


class FragmentSupportTests(unittest.TestCase):
    def setUp(self):
        self.a=record('r/1','r/2',1,'p',10,160,False)
        self.b=record('r/2','r/1',2,'p',210,360,True)

    def test_pair_is_one_fragment(self):
        result=list(collect([self.b,self.a],1000))
        self.assertEqual(len(result),1)
        self.assertEqual(result[0]['candidates']['p'][0]['insert_size'],350)

    def test_duplicate_alignment_does_not_increase_support(self):
        result=list(collect([self.a,self.a,self.b,self.b],1000))[0]
        self.assertEqual(len(result['candidates']['p']),1)

    def test_same_placement_keeps_only_best_score(self):
        worse=dict(self.a,score=100)
        result=list(collect([worse,self.a,self.b],1000))[0]
        self.assertEqual(len(result['candidates']['p']),1)
        self.assertEqual(result['candidates']['p'][0]['alignment_score'],320)

    def test_same_read_names_distinguished_by_mate_links(self):
        self.a.update(name='r',fragment_next={'name':'r'})
        self.b.update(name='r',fragment_prev={'name':'r'})
        self.assertEqual(list(collect([self.a,self.b],1000))[0]['status'],'supported')

    def test_different_paths_do_not_form_pair(self):
        self.b['placements'][0]['path']='q'
        self.assertEqual(list(collect([self.a,self.b],1000))[0]['status'],'no_concordant_path')

    def test_same_orientation_rejected(self):
        self.b['placements'][0]['reverse']=False
        self.assertFalse(list(collect([self.a,self.b],1000))[0]['candidates'])

    def test_missing_mate_is_retained(self):
        self.assertEqual(list(collect([self.a],1000))[0]['status'],'missing_mate')

    def test_long_insert_rejected(self):
        self.assertFalse(list(collect([self.a,self.b],300))[0]['candidates'])


if __name__ == '__main__':
    unittest.main()
