""" Unit tests for raster functions """

import unittest
import os
import sys
import numpy as np
from test_helper import md5sum, TESTDATA
import karta

from ..raster import _dem

class RegularGrid(unittest.TestCase):

    def setUp(self):
        pe = karta.raster.peaks(n=49)
        self.rast = karta.grid.RegularGrid(hdr={'nx':49, 'ny':49,
                                                'xllcorner':0.0,
                                                'yllcorner':0.0,
                                                'dx':30.0,
                                                'dy':30.0,
                                                'nbands':1}, Z=pe)
        return

    def test_add_rgrid(self):
        rast2 = karta.grid.RegularGrid(self.rast.get_hdr(),
                                       Z=np.random.random(self.rast.data.shape)) 
        res = self.rast + rast2
        self.assertTrue(np.all(res.data == self.rast.data+rast2.data))
        return

    def test_sub_rgrid(self):
        rast2 = karta.grid.RegularGrid(self.rast.get_hdr(),
                                       Z=np.random.random(self.rast.data.shape)) 
        res = self.rast - rast2
        self.assertTrue(np.all(res.data == self.rast.data-rast2.data))
        return

    def test_center_coords(self):
        ans = np.meshgrid(np.arange(15.0, 1471.0, 30.0),
                          np.arange(15.0, 1471.0, 30.0))
        self.assertEqual(0.0, np.sum(self.rast.center_coords()[0] - ans[0]))
        self.assertEqual(0.0, np.sum(self.rast.center_coords()[1] - ans[1]))
        return

    def test_resample(self):
        small = karta.raster.peaks(n=7)
        rast = self.rast.copy()
        rast.resample(210.0, 210.0)
        self.assertEqual(0.0, np.sum(rast.data - small))
        return

    def test_vertex_coords(self):
        ans = np.meshgrid(np.arange(0.0, 1471.0, 30.0),
                          np.arange(0.0, 1471.0, 30.0))
        self.assertEqual(0.0, np.sum(self.rast.vertex_coords()[0] - ans[0]))
        self.assertEqual(0.0, np.sum(self.rast.vertex_coords()[1] - ans[1]))
        return

    def test_get_region(self):
        creg = (15.0, 1485.0, 15.0, 1485.0)
        ereg = (0.0, 1500.0, 0.0, 1500.0)
        self.assertEqual(self.rast.get_region(reference='center'), creg)
        self.assertEqual(self.rast.get_region(reference='edge'), ereg)
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

    def test_aairead(self):
        grid = karta.grid.aairead(os.path.join(TESTDATA,'peaks49.asc'))
        self.assertTrue(False not in (grid.data == self.rast.data))
        return


class TestStructuredGrid(unittest.TestCase):

    def setUp(self):
        ii = np.arange(50.0)
        jj = np.arange(50.0)
        X, Y = np.meshgrid(np.sin(ii/25.0 * 2*np.pi),
                           np.sin(jj/50.0 * 2*np.pi))
        Z = karta.raster.witch_of_agnesi(50, 50)
        self.rast = karta.grid.StructuredGrid(X=X, Y=Y, Z=Z)

    def test_hdr(self):
        hdr = self.rast.get_hdr()
        self.assertEqual(hdr, {'xllcorner': 0.0,
                               'yllcorner': -0.12533323356430465,
                               'nbands':1})
        return

    def test_add_sgrid(self):
        rast2 = karta.grid.StructuredGrid(X=self.rast.X, Y=self.rast.Y,
                                          Z=np.random.random(self.rast.data.shape))
        res = self.rast + rast2
        self.assertTrue(np.all(res.data == self.rast.data+rast2.data))
        return

    def test_sub_sgrid(self):
        rast2 = karta.grid.StructuredGrid(X=self.rast.X, Y=self.rast.Y,
                                          Z=np.random.random(self.rast.data.shape))
        res = self.rast - rast2
        self.assertTrue(np.all(res.data == self.rast.data-rast2.data))
        return

class TestAAIGrid(unittest.TestCase):

    def setUp(self):
        pe = karta.raster.peaks(n=49)
        self.rast = karta.aaigrid.AAIGrid(pe, hdr={'ncols':49, 'nrows':49,
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

class TestDEMDriver(unittest.TestCase):

    def setUp(self):
        pass

    def test_nreps_re(self):
        nr = _dem.nreps_re("2(I4,I2,F7.4)")
        self.assertEqual(nr, 2)
        nr = _dem.nreps_re("2(3I6)")
        self.assertEqual(nr, 6)
        nr = _dem.nreps_re("4(6F3.7)")
        self.assertEqual(nr, 24)
        return

    def test_parse(self):
        p = _dem.parse("2(I4,I2,F7.4)", ' -81 0 0.0000  82 0 0.0000')
        self.assertEqual(p, [-81, 0, 0.0, 82, 0, 0.0])
        return

    def test_dtype(self):
        t = _dem.dtype("2(I4,D24.15,F7.4)")
        self.assertEqual(t, [int, _dem.coerce_float, _dem.coerce_float])
        return

    def test_reclen(self):
        n = _dem.reclen("2(I4,I2,F7.4)")
        self.assertEqual(n, [4, 2, 7])
        return

#class TestInterpolation(unittest.TestCase):
#
#    def test_idw(self):
#        # Test inverse weighted distance interpolation
#        Xm, Ym = np.meshgrid(np.arange(10), np.arange(10))
#        G = np.c_[Xm.flat, Ym.flat]
#        obs = np.array(((2,4,5),(6,3,4),(3,9,9)))
#        D = raster.interpolation.interp_idw(G, obs, wfunc=lambda d: 1.0/(d+1e-16)**0.5)
#        self.assertEqual(D, np.array([ 5.77099225,  5.73080888,  5.68934588,
#                            5.64656361,  5.60263926, 5.56391335,  5.54543646,
#                            5.55810434,  5.59468599,  5.64001991, 5.7666018 ,
#                            5.71712039,  5.66733266,  5.61528705,  5.55254697,
#                            5.48200769,  5.44194711,  5.4718629 ,  5.54176929,
#                            5.61255561, 5.76627041,  5.70140932,  5.64212525,
#                            5.59011147,  5.50963694, 5.36912109,  5.24490032,
#                            5.34965772,  5.49388958,  5.59812945, 5.7765376 ,
#                            5.67823283,  5.58875282,  5.57796102,  5.51352082,
#                            5.30130025,  4.00000002,  5.2696964 ,  5.49251053,
#                            5.61373077, 5.81944097,  5.68356449,  5.00000001,
#                            5.61034065,  5.60629104, 5.47740263,  5.34297441,
#                            5.43952794,  5.57499849,  5.6693856 , 5.91915228,
#                            5.83714501,  5.76171956,  5.78893204,  5.78328667,
#                            5.71759216,  5.65687878,  5.66124481,  5.70678887,
#                            5.75517987, 6.06116315,  6.0515271 ,  6.04980015,
#                            6.05331942,  6.02324779, 5.95580534,  5.88858853,
#                            5.8514621 ,  5.84418368,  5.85242989, 6.21638838,
#                            6.27774217,  6.35281161,  6.38711869,  6.31999906,
#                            6.19926611,  6.09004119,  6.01415787,  5.96971255,
#                            5.94708184, 6.35785785,  6.4954344 ,  6.70982294,
#                            6.88076712,  6.68090612, 6.43193841,  6.26024319,
#                            6.14772366,  6.07568133,  6.03068086, 6.45468497,
#                            6.63857194,  6.9936249 ,  8.99999996,  6.97005635,
#                            6.58706897,  6.37686191,  6.24425398,  6.15677403,
#                            6.09829941]))
#        return


if __name__ == "__main__":
    unittest.main()

