import unittest
import json
import tempfile
from pathlib import Path
from unittest.mock import patch
from advance_validation_graph import completed,advance


class DispatchTests(unittest.TestCase):
    def test_every_expected_task_must_succeed(self):
        for output in ('123_1|COMPLETED|0:0\n',
                       '123_1|COMPLETED|0:0\n123_2|FAILED|1:0\n',
                       '123_[1-2]|PENDING|0:0\n'):
            with patch('advance_validation_graph.remote',return_value=output):
                self.assertFalse(completed('123',1,2))
        with patch('advance_validation_graph.remote',return_value='123_2|COMPLETED|0:0\n123_1|COMPLETED|0:0\n'):
            self.assertTrue(completed('123',1,2))

    def test_projection_requires_matching_success_and_join_waits_for_all(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'work').mkdir()
            ledger=root/'validation-graph-launch.json'
            ledger.write_text(json.dumps({'freeze_sha256':'frozen',
                'jobs':{'cohort':{'calibration':'100','mapping':'123','native':'122'}}}))
            with patch('advance_validation_graph.HERE',root),patch(
                    'advance_validation_graph.completed',return_value=False),patch(
                    'advance_validation_graph.submit',return_value='124') as submit:
                advance()
                options=submit.call_args.args[-1]
                self.assertIn('--dependency=aftercorr:123',options)
                self.assertIn('--array=16-1343%300',options)
                self.assertEqual(json.loads(ledger.read_text())['jobs']['cohort']['evidence'],'124')
                submit.reset_mock()
                advance()
                submit.assert_not_called()
