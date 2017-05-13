import unittest
import os
import karta
import numpy as np
from test_helper import TESTDATA, TMPDATA

class AAITests(unittest.TestCase):

    def test_read_aai_corner(self):
        pe = peaks(n=49)
        control = karta.RegularGrid((0.0, 0.0, 30.0, 30.0, 0.0, 0.0), values=pe)
        grid = karta.read_aai(os.path.join(TESTDATA, "peaks49_corner.asc"))
        self.assertTrue(np.allclose(grid[::-1], control[:,:]))

    def test_read_aai_center(self):
        pe = peaks(n=49)
        control = karta.RegularGrid((0.0, 0.0, 30.0, 30.0, 0.0, 0.0), values=pe)
        grid = karta.read_aai(os.path.join(TESTDATA, "peaks49_center.asc"))
        self.assertTrue(np.allclose(grid[::-1], control[:,:]))

    def test_write_aai(self):
        pe = peaks(n=49)
        grid = karta.RegularGrid((0.0, 0.0, 30.0, 30.0, 0.0, 0.0), values=pe)
        ret = grid.to_aai(os.path.join(TMPDATA, "peaks49.asc"))
        self.assertEqual(ret, grid)

        with open(os.path.join(TMPDATA, "peaks49.asc")) as f:
            lines = f.readlines()
        self.assertEqual(len(lines), 55)
        self.assertEqual("".join(lines[:6]), """NCOLS 49
NROWS 49
XLLCORNER 0.0
YLLCORNER 0.0
CELLSIZE 30.0
NODATA_VALUE -9999
""")

    def test_write_aai_centered(self):
        pe = peaks(n=49)
        grid = karta.RegularGrid((0.0, 0.0, 30.0, 30.0, 0.0, 0.0), values=pe)
        grid.to_aai(os.path.join(TMPDATA, "peaks49.asc"), reference="center")
        with open(os.path.join(TMPDATA, "peaks49.asc")) as f:
            lines = f.readlines()
        self.assertEqual(len(lines), 55)
        self.assertEqual("".join(lines[:6]), """NCOLS 49
NROWS 49
XLLCENTER 15.0
YLLCENTER 15.0
CELLSIZE 30.0
NODATA_VALUE -9999
""")

def peaks(n=49):
    """ 2d peaks function of MATLAB logo fame. """
    X, Y = np.meshgrid(np.linspace(-3, 3, n), np.linspace(-3, 3, n))
    return 3.0 * (1-X)**2 * np.exp(-X**2 - (Y+1)**2) \
            - 10.0 * (X/5.0 - X**3 - Y**5) * np.exp(-X**2 - Y**2) \
            - 1.0/3.0 * np.exp(-(X+1)**2 - Y**2)

if __name__ == "__main__":
    unittest.main()
