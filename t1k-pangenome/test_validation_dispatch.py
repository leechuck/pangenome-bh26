import unittest
from unittest.mock import patch
from advance_validation_graph import completed


class DispatchTests(unittest.TestCase):
    def test_every_expected_task_must_succeed(self):
        for output in ('123_1|COMPLETED|0:0\n',
                       '123_1|COMPLETED|0:0\n123_2|FAILED|1:0\n',
                       '123_[1-2]|PENDING|0:0\n'):
            with patch('advance_validation_graph.remote',return_value=output):
                self.assertFalse(completed('123',1,2))
        with patch('advance_validation_graph.remote',return_value='123_2|COMPLETED|0:0\n123_1|COMPLETED|0:0\n'):
            self.assertTrue(completed('123',1,2))
