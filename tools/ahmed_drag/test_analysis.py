"""Analytic extrapolation checks; no restricted fixtures or CFD runs."""
import pathlib
import sys
import unittest

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import analyze


class AnalysisTest(unittest.TestCase):
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
