import unittest,random
from graft import graft,cds_of,rc

def random_seq(rng,n):return ''.join(rng.choice('ACGT') for _ in range(n))

class GraftTests(unittest.TestCase):
    def setUp(self):
        rng=random.Random(7)
        self.flank=random_seq(rng,50);self.e1=random_seq(rng,30);self.intr=random_seq(rng,40);self.e2=random_seq(rng,21);self.tail=random_seq(rng,50)
        self.seq=self.flank+self.e1+self.intr+self.e2+self.tail
        a=len(self.flank);b=a+30;c=b+40;d=c+21
        self.iv=[(a,b),(c,d)]
    def test_identity(self):
        s,iv=graft(self.seq,self.iv,'+',self.e1+self.e2)
        self.assertEqual(s,self.seq);self.assertEqual(iv,self.iv)
    def test_substitution_keeps_introns_and_flanks(self):
        new=self.e1[:10]+('A' if self.e1[10]!='A' else 'C')+self.e1[11:]+self.e2
        s,iv=graft(self.seq,self.iv,'+',new)
        self.assertEqual(cds_of(s,iv,'+'),new)
        self.assertEqual(s[:50],self.flank);self.assertEqual(s[iv[0][1]:iv[1][0]],self.intr);self.assertEqual(s[iv[1][1]:],self.tail)
    def test_indels_round_trip(self):
        new=self.e1[:5]+'GGG'+self.e1[5:]+self.e2[3:]
        s,iv=graft(self.seq,self.iv,'+',new)
        self.assertEqual(cds_of(s,iv,'+'),new);self.assertEqual(len(s),len(self.seq)+3-3)
        self.assertIn(self.intr,s)
    def test_minus_strand(self):
        seq=rc(self.seq);L=len(seq);iv=sorted((L-b,L-a) for a,b in self.iv)
        old=cds_of(seq,iv,'-');self.assertEqual(old,self.e1+self.e2)
        new=old[:12]+'T'+old[13:40]+'ACG'+old[40:]
        s,niv=graft(seq,iv,'-',new)
        self.assertEqual(cds_of(s,niv,'-'),new)
        self.assertTrue(s.endswith(rc(self.flank)))

if __name__=='__main__':unittest.main()
