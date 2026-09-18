import unittest
from build_candidates import norm2,t1k_backbones,partial_near,near

class CandidateTests(unittest.TestCase):
    def test_norm2(self):
        self.assertEqual(norm2('HLA-DRB1*08:03:02'),'08:03');self.assertEqual(norm2('A*01:01'),'01:01');self.assertIsNone(norm2('A*01'))
    def test_t1k_backbones_leave_one_out_and_label_match(self):
        labels={'X.1':{'A':dict(donor_id='X',cds_complete='1',exact_cds_alleles='A*02:01:01:01;A*02:01:01:02',exact_protein_alleles='')},
                'Y.1':{'A':dict(donor_id='Y',cds_complete='1',exact_cds_alleles='A*02:01:01:01',exact_protein_alleles='')},
                'Z.2':{'A':dict(donor_id='Z',cds_complete='0',exact_cds_alleles='',exact_protein_alleles='A*11:01:01')},
                'W.1':{'A':dict(donor_id='W',cds_complete='1',exact_cds_alleles='A*24:02:01',exact_protein_alleles='')}}
        hap={'X.1':'AAAA','Y.1':'AAAA','Z.2':'CCCC','W.1':'GGGG'}
        bb=t1k_backbones(['HLA-A*02:01:01','A*11:01'],'A',labels,donor='X',hap=hap)
        self.assertEqual(bb,['Y.1','Z.2'])   # X excluded (own donor), Y kept as distinct sequence, W label mismatch
    def test_partial_near_replaces_span(self):
        import random;rng=random.Random(3);cds='ATG'+''.join(rng.choice('ACGT') for _ in range(180))+'TAA'
        span=cds[40:100];mut=span[:10]+('A' if span[10]!='A' else 'C')+span[11:]
        rows=[dict(sequence=mut,alleles='G*01:02;G*01:02:01'),dict(sequence=span,alleles='G*01:01')]
        out=partial_near(cds,rows,3)
        self.assertEqual(len(out),1)
        new,(rep,d,alist)=next(iter(out.items()))
        self.assertEqual(rep,'G*01:02');self.assertEqual(new,cds[:40]+mut+cds[100:])
    def test_near_excludes_identical(self):
        cds='ATGAAACCCGGGTTTTAA'
        rows=[dict(sequence=cds,alleles='G*01:01'),dict(sequence=cds[:5]+'T'+cds[6:],alleles='G*01:02')]
        self.assertEqual([v[0] for v in near(cds,rows,3).values()],['G*01:02'])

if __name__=='__main__':unittest.main()
