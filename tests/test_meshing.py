import unittest

import numpy as np

from porphyry_gen.meshing import is_closed, marching_squares, marching_tetrahedra, mesh_volume


class MeshingTests(unittest.TestCase):
    def test_sphere_is_closed_with_correct_volume(self):
        g = np.arange(-50, 51, 5.0)
        X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
        v, t = marching_tetrahedra(40.0 - np.sqrt(X ** 2 + Y ** 2 + Z ** 2), (-50, -50, -50), (5, 5, 5))
        self.assertTrue(is_closed(t))
        self.assertAlmostEqual(mesh_volume(v, t) / (4 / 3 * np.pi * 40 ** 3), 1.0, delta=0.02)

    def test_body_touching_grid_edge_is_capped(self):
        g = np.arange(0, 101, 10.0)
        X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
        v, t = marching_tetrahedra(np.minimum(55.0 - Z, 1.0), (0, 0, 0), (10, 10, 10))
        self.assertTrue(is_closed(t))
        self.assertGreater(mesh_volume(v, t), 0)

    def test_circle_contour_closed(self):
        g = np.arange(-20, 21, 1.0)
        X, Y = np.meshgrid(g, g, indexing="ij")
        lines = marching_squares(10.0 - np.hypot(X, Y))
        self.assertEqual(len(lines), 1)
        pts, closed = lines[0]
        self.assertTrue(closed)
        r = np.hypot(pts[:, 0] - 20, pts[:, 1] - 20)
        self.assertTrue(np.allclose(r, 10.0, atol=0.3))


if __name__ == "__main__":
    unittest.main()
