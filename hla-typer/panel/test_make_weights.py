import unittest
from make_weights import tile
W=dict(exon=1,intron=0.1,intergenic=0.005)

class TileTests(unittest.TestCase):
    def test_tiles_whole_locus_without_gaps(self):
        rows=tile(100,[(20,60,[(20,30),(50,60)])],W)
        self.assertEqual(rows,[[0,20,0.005],[20,30,1],[30,50,0.1],[50,60,1],[60,100,0.005]])
    def test_two_genes_and_no_genes(self):
        self.assertEqual(tile(10,[],W),[[0,10,0.005]])
        rows=tile(100,[(10,20,[(10,12)]),(50,70,[(60,70)])],W)
        self.assertEqual(rows[0],[0,10,0.005]);self.assertEqual(rows[-1],[70,100,0.005])
        self.assertEqual(sum(b-a for a,b,_ in rows),100)
    def test_clipped_gene_stays_inside_locus(self):
        rows=tile(50,[(0,30,[(0,5)])],W)
        self.assertEqual(rows[0],[0,5,1])

if __name__=='__main__':unittest.main()
