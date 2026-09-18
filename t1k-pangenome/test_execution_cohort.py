import unittest
from score_linear_controls import verify_cohort_metadata


class ExecutionCohortTests(unittest.TestCase):
    def setUp(self):
        self.rows = [dict(donor='HG1',family='f1',fold='0',cram='/old/HG1.cram'),
                     dict(donor='HG2',family='f2',fold='1',cram='/old/HG2.cram')]

    def test_source_relocation_allowed(self):
        execution = [dict(row,cram='https://example.org/'+row['donor']+'.cram') for row in self.rows]
        verify_cohort_metadata(self.rows,execution)

    def test_selection_and_order_changes_rejected(self):
        for execution in ([dict(self.rows[0],family='changed'),self.rows[1]],
                          list(reversed(self.rows)),self.rows[:1],
                          [self.rows[0],self.rows[0]]):
            with self.assertRaises(ValueError):
                verify_cohort_metadata(self.rows,execution)

    def test_missing_metadata_rejected(self):
        row=dict(self.rows[0]);del row['fold']
        with self.assertRaises(ValueError):
            verify_cohort_metadata(self.rows,[row,self.rows[1]])
