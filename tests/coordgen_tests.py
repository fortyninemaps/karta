import unittest
import numpy as np
from karta.raster.coordgen import CoordinateGenerator
from karta.crs import CartesianCRS, LonLatWGS84, WebMercator

class CoordinateGeneratorTests(unittest.TestCase):

    def test_index_same_crs(self):
        cg = CoordinateGenerator([0, 0, 1, 1, 0, 0], (16, 16), CartesianCRS, CartesianCRS)
        x, y = cg[7, 5]
        self.assertEqual(x, 5.5)
        self.assertEqual(y, 7.5)
        return

    def test_slice_same_crs(self):
        cg = CoordinateGenerator([0, 0, 1, 1, 0, 0], (16, 16), CartesianCRS, CartesianCRS)
        x, y = cg[7:, 5:12:2]
        self.assertEqual(x.shape, (9 ,4))
        self.assertTrue(np.allclose(x[0,:], np.array([5.5, 7.5, 9.5, 11.5])))
        self.assertTrue(np.allclose(y[:,0], np.arange(7.5, 16, 1.0)))
        return

    def test_index_diff_crs(self):
        cg = CoordinateGenerator([0, 0, 1, 1, 0, 0], (16, 16), LonLatWGS84, WebMercator)
        x, y = cg[7, 5]
        self.assertAlmostEqual(x, 612257.1993630034)
        self.assertAlmostEqual(y, 837290.7323422873)
        return

    def test_slice_diff_crs(self):
        cg = CoordinateGenerator([0, 0, 1, 1, 0, 0], (16, 16), LonLatWGS84, WebMercator)
        x, y = cg[7:, 5:12:2]
        self.assertEqual(x.shape, (9 ,4))
        return

if __name__ == "__main__":
    unittest.main()
