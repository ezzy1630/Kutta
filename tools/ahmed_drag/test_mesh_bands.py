"""Analytic polygon checks; no generated meshes or restricted fixtures."""
import pathlib
import sys
import unittest

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import mesh_bands


class MeshBandsTest(unittest.TestCase):
    def test_cut_preserves_constant_pressure_load_and_exact_band_area(self):
        # A centroid rule assigns all 60 mm to a 40 mm band.
        face = np.array([[0., 0, 0], [.06, 0, 0], [.06, 1, 0], [0, 1, 0]])
        area, centre = mesh_bands.polygon_geometry(face)
        self.assertAlmostEqual(area, .06)
        np.testing.assert_allclose(centre, [.03, .5, 0])
        leading = mesh_bands.clipped_area(face, .04)
        remaining = mesh_bands.clipped_area(face, .04, keep_left=False)
        self.assertAlmostEqual(leading, .04)
        self.assertAlmostEqual(remaining, .02)
        cp_normal = .37
        self.assertAlmostEqual(cp_normal*(leading+remaining), cp_normal*area)

    def test_inside_outside_and_on_cut_vertices(self):
        triangle = np.array([[0., 0, 0], [1., 0, 0], [0, 1., 0]])
        self.assertAlmostEqual(mesh_bands.clipped_area(triangle, 2), .5)
        self.assertEqual(mesh_bands.clipped_area(triangle, -1), 0)
        self.assertEqual(mesh_bands.clipped_area(triangle, 0), 0)
        self.assertAlmostEqual(mesh_bands.clipped_area(triangle, .5), .375)


if __name__ == '__main__':
    unittest.main()
