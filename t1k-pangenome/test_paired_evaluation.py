import unittest
from paired_evaluation import evaluate,GENES,interval


def fixture():
    cohort=[dict(donor='d1',family='f1',stratum='EAS'),dict(donor='d2',family='f1',stratum='EAS'),dict(donor='d3',family='f2',stratum='AFR')]
    rows=[dict(method=m,donor=d['donor'],gene=g,fields=f,state='complete',eligible=1,called=1,
               correct=int(m=='new' or f==2))
          for m in ('T1K','new') for d in cohort for g in GENES for f in (2,4)]
    return rows,cohort


class PairedEvaluationTests(unittest.TestCase):
    def test_family_clustering_and_clear_gain(self):
        rows,cohort=fixture();r=evaluate(rows,cohort,'new',repetitions=100)
        self.assertTrue(r['protocol_pass'])
        self.assertEqual(r['contrasts']['ALL/4']['interval']['families'],2)
        self.assertEqual(r['contrasts']['ALL/4']['interval']['lower_95'],1)

    def test_pending_missing_duplicate_and_outside_donor_rejected(self):
        rows,cohort=fixture()
        for variant in (rows[:-1],rows+[rows[0]], [dict(rows[0],state='running')]+rows[1:],
                        [dict(rows[0],donor='foreign')]+rows[1:]):
            with self.assertRaises(ValueError):evaluate(variant,cohort,'new',repetitions=100)

    def test_ineligible_truth_and_failed_predictions(self):
        rows,cohort=fixture()
        for row in rows:
            if row['gene']=='A':row.update(eligible=0,called=0,correct=0)
            if row['method']=='new' and row['donor']=='d1':row.update(state='failed',called=0,correct=0)
        r=evaluate(rows,cohort,'new',repetitions=100)
        self.assertEqual(r['contrasts']['ALL/4']['interval']['eligible'],21)
        self.assertEqual(r['contrasts']['ALL/4']['failed_donors']['new'],1)
        self.assertFalse(r['protocol_pass'])

    def test_changed_eligibility_or_failed_credit_rejected(self):
        rows,cohort=fixture()
        rows[0].update(eligible=0,correct=0)
        with self.assertRaises(ValueError):evaluate(rows,cohort,'new',repetitions=100)
        rows,cohort=fixture();rows[-1]['state']='failed'
        with self.assertRaises(ValueError):evaluate(rows,cohort,'new',repetitions=100)

    def test_bootstrap_is_order_invariant(self):
        a={'b':(1,8),'a':(-1,4),'c':(3,8)}
        b=dict(reversed(list(a.items())))
        self.assertEqual(interval(a,100),interval(b,100))
