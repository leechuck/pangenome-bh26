import unittest
from validation_handoff import METHOD_DIRS,verify_grid


class HandoffTests(unittest.TestCase):
    def test_complete_grid_required_before_readiness(self):
        rows=[dict(donor=d,method=m,status='complete') for d in ('a','b') for m in METHOD_DIRS]
        self.assertTrue(verify_grid(rows,['a','b'])[1])
        rows[-1]['status']='running'
        self.assertFalse(verify_grid(rows,['a','b'])[1])
        rows[-1]['status']='failed'
        self.assertFalse(verify_grid(rows,['a','b'])[1])
        with self.assertRaises(ValueError):verify_grid(rows[:-1],['a','b'])
        rows[-1]=rows[0]
        with self.assertRaises(ValueError):verify_grid(rows,['a','b'])
