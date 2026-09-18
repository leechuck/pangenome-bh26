import tempfile
from pathlib import Path
import unittest
from graph_recruitment import select_reads,paired


def write(path,ids,base='A'):
    path.write_text(''.join('@'+i+'\n'+base*4+'\n+\nIIII\n' for i in ids))


class RecruitmentTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.raw=[self.root/'raw1.fq',self.root/'raw2.fq']
        self.native=[self.root/'native1.fq',self.root/'native2.fq']
        for i,(raw,native) in enumerate(zip(self.raw,self.native),1):
            write(raw,[f'{r}/{i}' for r in ('a','b','c')]);write(native,[f'a/{i}'])

    def test_union_deduplicates_paths_and_preserves_native(self):
        result=select_reads(self.raw,self.native,self.root,['b','b','a'])
        self.assertEqual(result['output_pairs'],2);self.assertEqual(result['added_pairs'],1)
        self.assertEqual([r[0] for r in paired(self.root/'r1.fq',self.root/'r2.fq')],['a','b'])

    def test_all_control_and_native_parity_selection(self):
        result=select_reads(self.raw,self.native,self.root)
        self.assertEqual(result['output_pairs'],1)
        result=select_reads(self.raw,self.native,self.root,all_reads=True)
        self.assertEqual(result['output_pairs'],3)

    def test_foreign_graph_fragment_rejected(self):
        with self.assertRaisesRegex(ValueError,'absent'):
            select_reads(self.raw,self.native,self.root,['other'])

    def test_changed_native_read_rejected(self):
        write(self.native[0],['a/1'],base='T')
        with self.assertRaisesRegex(ValueError,'differs'):
            select_reads(self.raw,self.native,self.root)

    def test_unpaired_and_duplicate_reads_rejected(self):
        write(self.raw[1],['a/2','b/2'])
        with self.assertRaisesRegex(ValueError,'Unequal'):
            select_reads(self.raw,self.native,self.root)
        for i,p in enumerate(self.raw,1):write(p,[f'a/{i}',f'a/{i}'])
        with self.assertRaisesRegex(ValueError,'Duplicate raw'):
            select_reads(self.raw,self.native,self.root)
