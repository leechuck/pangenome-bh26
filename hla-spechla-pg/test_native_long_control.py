import unittest
from native_long_control import patch_insertion_phaser

class InsertionCompatibilityTests(unittest.TestCase):
    def test_preserves_inputs_and_inference_settings(self):
        s='before\nSpecHap --ncs --window_size 15000 -N --vcf input.vcf.gz --frag reads --out result\nafter\n'
        t=patch_insertion_phaser(s)
        self.assertEqual(t.replace('--protocols nanopore --weights 1','-N'),s)
    def test_refuses_unknown_or_duplicate_source(self):
        for s in ('changed source',('SpecHap --ncs --window_size 15000 -N --vcf x\n')*2):
            with self.assertRaises(ValueError):patch_insertion_phaser(s)

class EmptyAssemblyTests(unittest.TestCase):
    def test_empty_contig_list_uses_no_contig_branch(self):
        import subprocess,tempfile
        from pathlib import Path
        from native_long_control import patch_empty_assembly
        source='if [ ! -f "$outdir/extract.fa" ]; then echo noreads; fi\nif [ ! -f "$outdir/id.list" ]; then echo no_contigs; else echo align; fi\n'
        # Production guards have the test and then on separate lines.
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'extract.fa').write_text('>r\nACGT\n');(p/'id.list').touch()
            patched=patch_empty_assembly(source)
            cmd='outdir='+d+'\n'+patched
            self.assertEqual(subprocess.check_output(['bash','-c',cmd],text=True).strip(),'no_contigs')
            (p/'id.list').write_text('contig1\n')
            self.assertEqual(subprocess.check_output(['bash','-c',cmd],text=True).strip(),'align')

if __name__=='__main__':unittest.main()
