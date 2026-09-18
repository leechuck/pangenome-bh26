import unittest
from map_personalized import sampling_decision


class SamplingTests(unittest.TestCase):
    def test_low_coverage_retains_unsampled_graph(self):
        self.assertEqual(sampling_decision({'genes':{'A':dict(usable=False,fallback_reason='low coverage')}},'A'),(None,'low coverage'))

    def test_selection_disabled_is_explicit_control(self):
        self.assertEqual(sampling_decision({},'A',True),(None,'selection_disabled_control'))

    def test_only_valid_measured_coverage_is_used(self):
        self.assertEqual(sampling_decision({'genes':{'A':dict(usable=True,personalization_coverage=69)}},'A'),(69,None))
        with self.assertRaises(ValueError):
            sampling_decision({'genes':{'A':dict(usable=True,personalization_coverage=None)}},'A')


if __name__ == '__main__':
    unittest.main()
