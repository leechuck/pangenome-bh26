import unittest
from panel_comparison import audit_membership, paired


class PanelTests(unittest.TestCase):
    def test_additive_rejects_missing_hprc_and_excluded_donors(self):
        meta = {'h': {'cohort': 'HPRC-Rest'}, 'a': {'cohort': 'APR'}}
        panels = dict(hprc=['h'], asian_matched=['a'], full=['h', 'a'])
        self.assertEqual(audit_membership(panels, meta, set())['added_haplotypes'], 1)
        with self.assertRaises(ValueError):
            audit_membership(panels, meta, {'a'})
        with self.assertRaises(ValueError):
            audit_membership(dict(panels, full=['a']), meta, set())

    def test_paired_counts_gains_and_losses_not_unpaired_totals(self):
        rows = [dict(method=m, donor=d, gene='A', level='two_field', stratum='EAS',
                     state='complete', eligible='1', correct=c)
                for m, d, c in [('a', 'd1', '1'), ('a', 'd2', '0'),
                                ('b', 'd1', '0'), ('b', 'd2', '1')]]
        r = paired(rows, 'a', 'b', 'two_field')
        self.assertEqual((r['gains'], r['losses'], r['eligible']), (1, 1, 2))
        rows[0]['state'] = 'running'
        self.assertIsNone(paired(rows, 'a', 'b', 'two_field'))


if __name__ == '__main__':
    unittest.main()
