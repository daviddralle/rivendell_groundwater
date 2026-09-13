"""Regression checks for the scientific accounting errors corrected in v0.2.0."""
import unittest
import numpy as np
import pandas as pd
from analysis_utils import interpolate_short_gaps, integrate_profile, surveyed_elevation
from neutron_probe_processing import assign_calibration_periods


class AccountingTests(unittest.TestCase):
    def test_only_short_bracketed_gaps_are_filled(self):
        dates = pd.date_range('2025-01-01', periods=15, freq='D')
        data = pd.DataFrame({'head':[np.nan, 1, np.nan, 3, np.nan, np.nan, np.nan,
                                   np.nan, 8, 9, np.nan, np.nan, 12, np.nan, np.nan]}, index=dates)
        actual = interpolate_short_gaps(data, 2)['head']
        self.assertEqual(actual.iloc[2], 2.)
        np.testing.assert_allclose(actual.iloc[10:12], [10., 11.])
        self.assertTrue(actual.iloc[4:8].isna().all())
        self.assertTrue(actual.iloc[[0,13,14]].isna().all())
        pd.testing.assert_frame_equal(data, pd.DataFrame({'head':[np.nan, 1, np.nan, 3, np.nan, np.nan, np.nan,
                                   np.nan, 8, 9, np.nan, np.nan, 12, np.nan, np.nan]}, index=dates))

    def test_integral_uses_fixed_bounds_and_millimeters(self):
        z = np.arange(0., 3.25, 0.25)
        wet, dry = 0.2 + 0.02*z, 0.1 + 0.01*z
        actual = integrate_profile(z, wet, .75, 2.5) - integrate_profile(z, dry, .75, 2.5)
        expected = (0.1*(2.5-.75)+.005*(2.5**2-.75**2))*1000
        self.assertAlmostEqual(actual, expected, places=10)
        # Extending a log changes neither the compared volume nor its storage.
        self.assertAlmostEqual(integrate_profile(z,wet,.75,2.5),
                               integrate_profile(z[:-1],wet[:-1],.75,2.5))
        with self.assertRaises(ValueError):
            integrate_profile(z[:9],wet[:9],.75,2.5)
        with self.assertRaises(ValueError):
            integrate_profile([0,.5,1.5,2],[.1,.1,.1,.1],.75,1.75)

    def test_calibration_matches_serial_and_period(self):
        observations = pd.DataFrame({'Date':pd.to_datetime(['2020-03-01','2020-03-01','2020-07-01','2020-03-01']),
                                     'Serial Number':['A','B','A','unknown']})
        intervals = pd.DataFrame({'Start':pd.to_datetime(['2020-01-01','2020-02-01','2020-07-01']),
                                 'End':pd.to_datetime(['2020-07-01','2021-01-01','2021-01-01']),
                                 'Serial Number':['A','B','A'], 'Conversion ID':['A1','B1','A2']})
        actual, rejected = assign_calibration_periods(observations, intervals)
        self.assertEqual(actual['Conversion ID'].tolist(), ['A1','B1','A2'])
        self.assertEqual(rejected['Serial Number'].tolist(), ['unknown'])

    def test_survey_datum_is_read_explicitly(self):
        self.assertEqual(surveyed_elevation('<td>Z</td>\\n<td>468.07716</td>'), 468.07716)
        with self.assertRaises(ValueError):
            surveyed_elevation('<td>X</td><td>468.07716</td>')


if __name__ == '__main__':
    unittest.main()
