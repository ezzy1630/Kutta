"""Analytic reporting regressions; no restricted data or CFD solves."""
import math
import pathlib
import sys
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import report


class AhmedReportTest(unittest.TestCase):
    def test_pressure_regions_exclude_side_taps(self):
        angle = math.radians(35)
        x = -.12
        z = .338 - (x + .222 * math.cos(angle)) * math.tan(angle)
        xyz = np.array([[x, 0, z], [x, -.1945, z - .02], [0, 0, .1]])
        walls = np.array([tuple(p) + (-.3,) for p in xyz],
                         dtype=[(k, 'f8') for k in ['x', 'y', 'z', 'cp']])
        with patch.object(report, 'lienhart', return_value=(xyz, np.full(3, -.3))):
            comparisons = report.cp_compare(walls, 35)
        for *_, regions in comparisons:
            np.testing.assert_array_equal(regions['slant'], [True, False, False])
            np.testing.assert_array_equal(regions['base'], [False, False, True])

    def test_25_degree_slant_and_requested_area(self):
        angle = math.radians(25)
        walls = np.array([(-.19, -math.sin(angle), -math.cos(angle), 1, -.5)],
                         dtype=[(x, 'f8') for x in ['x', 'nx', 'nz', 'area', 'cp']])
        with patch.object(report, 'walls', return_value=walls):
            parts, _ = report.ahmed_breakdown(
                {'dir': pathlib.Path('.'), 'point': '25', 'command': ['--area', '2']})
        self.assertAlmostEqual(parts['slant'], .5 * math.sin(angle) / 2)
        self.assertEqual(parts['front'], 0)
        self.assertEqual(parts['base'], 0)


if __name__ == '__main__':
    unittest.main()
