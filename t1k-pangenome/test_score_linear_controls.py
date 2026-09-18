import unittest
from score_linear_controls import prediction, compare


class ControlScoringTests(unittest.TestCase):
    def calls(self, alleles, unresolved=False, fields=4):
        return {'A':{'resolutions':{str(fields):[
            dict(alleles=alleles,unresolved=unresolved)]*2}}}

    def test_two_field_uses_existing_numeric_comparison_keys(self):
        pair = prediction(self.calls(['A*02:01'],fields=2),'A',2)
        self.assertEqual(compare(pair,[{'02:01'},{'02:01'}])[1],1)

    def test_wrong_alternative_is_not_rewarded(self):
        pair = prediction(self.calls(['A*02:01:01:01','A*02:02:01:01']),'A',4)
        self.assertEqual(compare(pair,[{'A*02:01:01:01'}]*2)[1],0)

    def test_unresolved_member_cannot_be_ignored(self):
        pair = prediction(self.calls(['A*02:01:01:01'],True),'A',4)
        self.assertEqual(compare(pair,[{'A*02:01:01:01'}]*2),(0,0,0))


if __name__ == '__main__':
    unittest.main()
