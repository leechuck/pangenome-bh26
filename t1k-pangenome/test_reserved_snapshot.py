import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from evidence_io import sha
from score_reserved_snapshot import run


class SnapshotTests(unittest.TestCase):
    def test_unready_snapshot_cannot_load_truth_or_create_unblinding_marker(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);validation=root/'validation';validation.mkdir()
            frozen=validation/'FROZEN_VALIDATION.json';frozen.write_text('{}')
            evaluation=validation/'EVALUATION_FREEZE.json'
            evaluation.write_text(json.dumps(dict(code_sha256={},input_sha256={},inference_freeze_sha256=sha(frozen))))
            (validation/'HANDOFF_FREEZE.json').write_text(json.dumps(dict(code_sha256={},evaluation_freeze_sha256=sha(evaluation))))
            snapshot=root/'snapshot';snapshot.mkdir()
            (snapshot/'EXECUTION_READY.json').write_text(json.dumps(dict(ready=False)))
            with patch('score_reserved_snapshot.HERE',root),patch('score_reserved_snapshot.truth_a',side_effect=AssertionError('Truth must remain unopened')):
                with self.assertRaisesRegex(ValueError,'not ready'):
                    run(snapshot,root/'output')
            self.assertFalse((validation/'UNBLINDED.json').exists())
