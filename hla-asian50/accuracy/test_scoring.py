import unittest,tempfile
from pathlib import Path
from score_panels import gt,predictions
class ScoringTests(unittest.TestCase):
    def test_sequence_identity_and_phase(self):
        self.assertEqual(gt('0|1',['A','AT']),gt('2/0',['AT','G','A']))
        self.assertEqual(gt('1/2',['A','AT','AT']),('AT','AT'))
        self.assertNotEqual(gt('0/1',['A','AT']),gt('0/1',['A','AG']))
    def test_missing_never_reference(self):
        for s in ['.','./.','0/.','.|1','0']:
            self.assertIsNone(gt(s,['A','T']))
    def test_vcf_format_and_missing(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x.vcf'
            p.write_text('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nctg1\t7\t.\tA\tT\t.\tPASS\t.\tGQ:GT\t90:1/0\nctg1\t9\t.\tC\tG\t.\tPASS\t.\tGT:GQ\t.:0\n')
            r=predictions(p);self.assertEqual(r[7,'A'][0],('A','T'));self.assertIsNone(r[9,'C'][0])
    def test_end_to_end_coverage_and_alt_remapping(self):
        import subprocess,sys,gzip,csv
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'source').mkdir();(root/'graph').mkdir()
            (root/'source/donors.tsv').write_text('donor\tstratum\tpopulation\tfold\nD\tEAS\tCHB\t0\n')
            header='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tD\n'
            ins='C'+'G'*50
            truth=header+f'ctg1\t100000\tx\tA\tT,C\t.\tPASS\t.\tGT\t0|1\nctg1\t100010\ty\tC\t{ins}\t.\tPASS\t.\tGT\t0|1\nctg1\t100100\tz\tG\tA\t.\tPASS\t.\tGT\t0|1\n'
            with gzip.open(root/'graph/full.top.vcf.gz','wt') as f:f.write(truth)
            full=header+f'ctg1\t100000\t.\tA\tC,T\t.\tPASS\t.\tGT\t2/0\nctg1\t100010\t.\tC\t{ins}\t.\tPASS\t.\tGT\t.\n'
            hprc=truth.replace('0|1','0/1')
            for arm,content in [('full',full),('hprc',hprc)]:
                out=root/f'pangenie/D/{arm}';out.mkdir(parents=True);(out/'calls_genotyping.vcf').write_text(content)
            subprocess.run([sys.executable,str(Path(__file__).with_name('score_panels.py')),'--run',str(root),'--out',str(root/'out')],check=True,capture_output=True)
            with open(root/'out/panel_per_donor.tsv') as f:rows=list(csv.DictReader(f,delimiter='\t'))
            def row(u,k,a):return next(r for r in rows if (r['universe'],r['variant_class'],r['arm'])==(u,k,a))
            self.assertEqual(row('shared_sites','SNV','full')['correct'],'1')
            self.assertEqual(row('all_truth','SNV','full')['n'],'2')
            self.assertEqual(row('all_truth','SNV','full')['site_present'],'1')
            self.assertEqual(row('shared_sites','truth_SV_length','full')['correct'],'0')
            self.assertEqual(row('shared_sites','truth_SV_length','hprc')['correct'],'1')
if __name__=='__main__':unittest.main()
