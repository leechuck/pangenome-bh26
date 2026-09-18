import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from evidence_io import sha
from run_hgsvc_stage import verify
from run_hgsvc_graph import coordinates, GENES


class ExecutionTest(unittest.TestCase):
    def test_tampered_code_and_reference_block_dispatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            (root/'code.py').write_text('original')
            (root/'reference.json').write_text('reference')
            frozen=root/'FROZEN_VALIDATION.json'
            frozen.write_text(json.dumps(dict(code_sha256={'code.py':sha(root/'code.py')},metadata_sha256={},asset_sha256={'reference.json':sha(root/'reference.json')})))
            digest=sha(frozen)
            with patch('run_hgsvc_stage.ROOT',root):
                verify(root,digest)
                (root/'code.py').write_text('changed')
                with self.assertRaises(ValueError):verify(root,digest)
                (root/'code.py').write_text('original')
                (root/'reference.json').write_text('changed')
                with self.assertRaises(ValueError):verify(root,digest)

    def test_complete_two_panel_28_donor_mapping(self):
        expected={(d,p,g) for d in range(28) for p in ('hprc','hprc_asian') for g in GENES}
        actual=[coordinates('mapping',i,28) for i in range(448)]
        self.assertEqual(set(actual),expected)
        self.assertEqual(len(actual),len(set(actual)))
        with self.assertRaises(ValueError):coordinates('mapping',448,28)


if __name__=='__main__':unittest.main()
