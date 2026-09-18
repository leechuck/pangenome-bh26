import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import advance_hgsvc_validation as dispatcher


class DispatchTest(unittest.TestCase):
    def test_mapping_requires_every_calibration_task(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            ledger=out/'EXECUTION_LAUNCH.json'
            ledger.write_text(json.dumps(dict(jobs={'calibration':'123'},freeze_sha256='digest',pending_submission=None)))
            complete='\n'.join('123_'+str(i)+'|COMPLETED|0:0' for i in range(28))
            with patch.object(dispatcher,'OUT',out),patch.object(dispatcher,'remote',return_value=complete.replace('123_27|COMPLETED','123_27|RUNNING')),patch.object(dispatcher,'submit') as submit:
                dispatcher.main()
                submit.assert_not_called()
            with patch.object(dispatcher,'OUT',out),patch.object(dispatcher,'remote',return_value=complete),patch.object(dispatcher,'submit',return_value='456') as submit:
                dispatcher.main()
                self.assertEqual(submit.call_count,1)
            self.assertEqual(json.loads(ledger.read_text())['jobs']['mapping'],'456')

    def test_uncertain_submission_is_not_retried(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            (out/'EXECUTION_LAUNCH.json').write_text(json.dumps(dict(pending_submission={'stage':'mapping'})))
            with patch.object(dispatcher,'OUT',out),patch.object(dispatcher,'submit') as submit:
                with self.assertRaises(RuntimeError):dispatcher.main()
                submit.assert_not_called()


if __name__=='__main__':unittest.main()
