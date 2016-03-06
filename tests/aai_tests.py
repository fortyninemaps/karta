import unittest
import karta
import numpy as np

class AAIGridTests(unittest.TestCase):

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

    def test_resize(self):
        orig = self.rast.data.copy()
        x0, x1, y0, y1 = self.rast.get_region()
        self.rast.resize((x0, x1/2.0, y0, y1/2.0))
        self.assertTrue(np.all(self.rast.data == orig[:25,:25]))
        return

if __name__ == "__main__":
    unittest.main()
