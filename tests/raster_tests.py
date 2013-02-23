""" Unit tests for raster functions """

import unittest
import os
import numpy as np
from test_helper import md5sum
import karta.raster as raster

class TestAAIGrid(unittest.TestCase):

    def setUp(self):
        pe = raster.raster.peaks(n=49)
        self.rast = raster.aaigrid.AAIGrid(pe, hdr={'ncols':49, 'nrows':49,
                                                    'xllcorner':0.0,
                                                    'yllcorner':0.0,
                                                    'cellsize':30.0,
                                                    'nodata_value':-9999})
        return

    def test_region_centered(self):
        reg = self.rast.get_region()
        self.assertEqual(reg, (15.0, 1485.0, 15.0, 1485.0))
        return

    def test_minmax(self):
        minmax = self.rast.minmax()
        self.assertEqual(minmax, (-6.5466445243204294, 8.075173545159231))
        return

    def test_get_indices(self):
        ind = self.rast.get_indices(0.0, 0.0)
        self.assertEqual(ind, (0, 48))
        ind = self.rast.get_indices(1485.0, 1485.0)
        self.assertEqual(ind, (48, 0))
        return

    def test_profile(self):
        prof = self.rast.get_profile([(0.0, 0.0), (1485.0, 1485.0)],
                                     resolution=30.0)
        expected = np.array([3.22353596e-05, 3.22353596e-05, 1.12058629e-04,
                             3.61838717e-04, 3.61838717e-04, 1.08395539e-03,
                             3.00817013e-03, 3.00817013e-03, 7.72015842e-03,
                             1.82832037e-02, 3.98497372e-02, 3.98497372e-02,
                             7.96679278e-02, 1.45455904e-01, 1.45455904e-01,
                             2.41123457e-01, 3.60000957e-01, 4.78441147e-01,
                             4.78441147e-01, 5.55779901e-01, 5.47061156e-01,
                             5.47061156e-01, 4.29133262e-01, 2.28899450e-01,
                             2.28899450e-01, 3.21030706e-02, -4.58152669e-02,
                             7.74741386e-02, 7.74741386e-02, 3.89805020e-01,
                             7.72078470e-01, 7.72078470e-01, 1.05657945e+00,
                             1.12541675e+00, 9.81011843e-01, 9.81011843e-01,
                             7.35348921e-01, 5.25639991e-01, 5.25639991e-01,
                             4.21139020e-01, 3.89628162e-01, 3.89628162e-01,
                             3.40764035e-01, 2.03711858e-01, -2.20472642e-02,
                             -2.20472642e-02, -2.72916549e-01, -4.67917982e-01,
                             -4.67917982e-01, -5.56069225e-01, -5.35273179e-01,
                             -4.40491855e-01, -4.40491855e-01, -3.18061354e-01,
                             -2.04492912e-01, -2.04492912e-01, -1.18122152e-01,
                             -6.16397937e-02, -2.91470271e-02, -2.91470271e-02,
                             -1.25010468e-02, -4.85698059e-03, -4.85698059e-03,
                             -1.70224503e-03, -5.33337096e-04, -5.33337096e-04,
                             -1.46658120e-04, -3.39563401e-05, -5.86418787e-06,
                             -5.86418787e-06, -5.86418787e-06])
        self.assertTrue(np.allclose(prof, expected))
        return

    def test_resize(self):
        orig = self.rast.data.copy()
        x0, x1, y0, y1 = self.rast.get_region()
        self.rast.resize((x0, x1/2.0, y0, y1/2.0))
        self.assertTrue(False not in (self.rast.data == orig[:25,:25]))
        return


if __name__ == "__main__":
    unittest.main()

