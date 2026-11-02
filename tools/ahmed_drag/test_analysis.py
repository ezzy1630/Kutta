"""Analytic extrapolation checks; no restricted fixtures or CFD runs."""
import pathlib
import sys
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import analyze


class AnalysisTest(unittest.TestCase):
    def test_pressure_errors_use_common_taps(self):
        rows = [{'run': name, 'mesh': {'slant_deg': 35.0}, 'cp': {}} for name in ('coarse', 'fine')]
        def comparison(path, angle):
            self.assertEqual(angle, 35)
            reference = np.zeros(3)
            kutta = np.array([100., 2., 3.]) if path == 'coarse' else np.array([100., 4., 5.])
            distance = np.array([.01, 0, 0]) if path == 'coarse' else np.array([0, 0, .01])
            groups = {region: np.ones(3, dtype=bool) for region in ('slant', 'base')}
            return [(name, None, reference, kutta, distance, groups) for name in ('experiment', 'CFX-5 SST')]
        with patch.object(analyze.report, 'walls', side_effect=lambda p: p), \
             patch.object(analyze.report, 'cp_compare', side_effect=comparison):
            analyze.common_pressure_comparisons(rows, ['coarse', 'fine'])
        self.assertEqual(rows[0]['cp_common']['experiment slant']['points'], 1)
        self.assertEqual(rows[1]['cp_common']['experiment slant']['points'], 1)
        self.assertEqual(rows[0]['cp_common']['experiment slant']['rms'], 2)
        self.assertEqual(rows[1]['cp_common']['experiment slant']['rms'], 4)

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
