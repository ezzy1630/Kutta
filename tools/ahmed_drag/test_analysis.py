"""Analytic extrapolation checks; no restricted fixtures or CFD runs."""
import pathlib
import json
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import analyze


class AnalysisTest(unittest.TestCase):
    def test_upstream_turbulence_is_volume_weighted_and_empty_is_unknown(self):
        cells = np.array([(-1.3, 0, .4, 1, .9, .01, .00015, 4, 1e-6),
                          (-1.4, 0, .5, 3, 1.1, .02, .0006, 8, 2e-6),
                          (-.5, 0, .4, 100, 50, 50, 1, 50, 50)],
                         dtype=[(name, 'f8') for name in
                                ('x', 'y', 'z', 'volume', 'u', 'p', 'k', 'omega', 'mu_t')])
        row = analyze.upstream_turbulence(cells, 1e-6)
        self.assertEqual(row['cells'], 2)
        self.assertAlmostEqual(row['tu_over_U_inlet'], .0175)
        self.assertAlmostEqual(row['mu_t_over_mu'], 1.75)
        self.assertAlmostEqual(row['u'], 1.05)
        self.assertAlmostEqual(row['cp_outlet_gauge'], .035)
        empty = analyze.upstream_turbulence(cells[2:], 1e-6)
        self.assertEqual(empty['cells'], 0)
        self.assertIsNone(empty['tu_over_U_inlet'])

    def test_rejects_unconfirmed_configuration_before_reading_fields(self):
        with tempfile.TemporaryDirectory() as directory:
            path = pathlib.Path(directory)
            (path / 'result.json').write_text(json.dumps({'configuration_verified': False}))
            with self.assertRaises(ValueError):
                analyze.summarize(path)

    def test_pressure_errors_use_common_taps(self):
        rows = [{'run': name, 'mesh': {'slant_deg': 35.0}, 'cp': {}} for name in ('coarse', 'fine')]
        def comparison(path, angle):
            self.assertEqual(angle, 35)
            reference = np.zeros(3)
            kutta = np.array([100., 2., 3.]) if path == 'coarse' else np.array([100., 4., 5.])
            distance = np.array([.01, 0, 0]) if path == 'coarse' else np.array([0, 0, .01])
            groups = {region: np.ones(3, dtype=bool) for region in ('slant', 'base')}
            start = -.222 * np.cos(np.radians(35))
            xyz = np.array([[start + .01, 0, .33], [start + .02, 0, .32], [-.01, 0, .22]])
            return [(name, xyz, reference, kutta, distance, groups) for name in ('experiment', 'CFX-5 SST')]
        with patch.object(analyze.report, 'walls', side_effect=lambda p: p), \
             patch.object(analyze.report, 'cp_compare', side_effect=comparison):
            analyze.common_pressure_comparisons(rows, ['coarse', 'fine'])
        self.assertEqual(rows[0]['cp_common']['experiment slant']['points'], 1)
        self.assertEqual(rows[1]['cp_common']['experiment slant']['points'], 1)
        self.assertEqual(rows[0]['cp_common']['experiment slant']['rms'], 2)
        self.assertEqual(rows[1]['cp_common']['experiment slant']['rms'], 4)
        self.assertEqual(rows[0]['cp_common']['experiment slant_leading_40mm']['points'], 1)
        self.assertEqual(rows[1]['cp_common']['experiment slant_leading_40mm']['kutta_mean'], 4)
        self.assertEqual(rows[0]['cp_common']['experiment slant_remaining']['points'], 0)
        self.assertIsNone(rows[0]['cp_common']['experiment slant_remaining']['rms'])

    def test_unequal_refinement_second_order(self):
        h = np.array([.04, .028, .0235])
        fit = analyze.richardson(h, .27 + 20 * h**2)
        self.assertTrue(fit['accepted'])
        self.assertAlmostEqual(fit['order'], 2)
        self.assertAlmostEqual(fit['extrapolated_cd'], .27)
        self.assertFalse(fit['asymptotic_verified'])

    def test_rejects_nonmonotonic_and_no_positive_order(self):
        self.assertFalse(analyze.richardson([4, 2, 1], [.3, .4, .3])['accepted'])
        self.assertFalse(analyze.richardson([4, 2, 1], [.4, .3, .1])['accepted'])


if __name__ == '__main__':
    unittest.main()
