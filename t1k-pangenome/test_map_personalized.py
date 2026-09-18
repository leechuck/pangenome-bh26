import unittest
from map_personalized import sampling_decision, fragment_parameters


class SamplingTests(unittest.TestCase):
    def test_calibration_requires_same_reads_and_valid_distribution(self):
        record = dict(status='complete', reads_sha256={'r1':'a','r2':'b'},
                      estimate=dict(usable=True,fragment_mean=400,fragment_stdev=60))
        self.assertEqual(fragment_parameters(record, record['reads_sha256']),
                         ['--fragment-mean','400','--fragment-stdev','60'])
        with self.assertRaises(ValueError):
            fragment_parameters(record, {'r1':'changed','r2':'b'})
        for invalid in (0, -1, float('nan'), float('inf'), True, '60'):
            record['estimate']['fragment_stdev'] = invalid
            with self.assertRaises(ValueError):
                fragment_parameters(record, record['reads_sha256'])

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
