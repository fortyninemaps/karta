import unittest
import os.path
import numpy as np
from test_helper import TMPDATA

import karta
from karta.raster import _gdal

class GdalTests(unittest.TestCase):

    def test_numpy_type_coercion(self):
        self.assertEqual(_gdal.numpy_dtype(2), np.uint16)
        self.assertEqual(_gdal.numpy_dtype(3), np.int16)
        self.assertEqual(_gdal.numpy_dtype(4), np.uint32)
        self.assertEqual(_gdal.numpy_dtype(5), np.int32)
        self.assertEqual(_gdal.numpy_dtype(6), np.float32)
        self.assertEqual(_gdal.numpy_dtype(7), np.float64)
        self.assertEqual(_gdal.numpy_dtype(8), np.complex64)
        self.assertEqual(_gdal.numpy_dtype(9), np.complex64)
        self.assertEqual(_gdal.numpy_dtype(10), np.complex64)
        self.assertEqual(_gdal.numpy_dtype(11), np.complex64)
        return

    def test_write_read(self):
        # try writing a file, then read it back in and verify that it matches
        v = peaks(500)[:100,:]
        utm7 = karta.crs.ProjectedCRS("+proj=utm +zone=7 +north +datum=WGS84",
                                      "UTM 7N (WGS 84)")
        g = karta.RegularGrid([15.0, 15.0, 30.0, 30.0, 0.0, 0.0], v, crs=utm7)

        fpath = os.path.join(TMPDATA, "test.tif")
        ret = g.to_geotiff(fpath, compress=None)
        self.assertEqual(ret, g)
        gnew = karta.read_geotiff(fpath)

        self.assertTrue("+proj=utm" in gnew.crs.get_proj4())
        self.assertTrue("+zone=7" in gnew.crs.get_proj4())
        self.assertEqual(g.transform, gnew.transform)
        self.assertEqual(g.values.dtype, gnew.values.dtype)
        self.assertTrue(np.all(g[:,:] == gnew[:,:]))
        return

    def test_write_read_disk(self):
        # try writing a file, then open it without loading into memory and verify
        v = peaks(500)[:100,:]
        utm7 = karta.crs.ProjectedCRS("+proj=utm +zone=7 +north +datum=WGS84",
                                      "UTM 7N (WGS 84)")
        g = karta.RegularGrid([15.0, 15.0, 30.0, 30.0, 0.0, 0.0], v, crs=utm7)

        fpath = os.path.join(TMPDATA, "test.tif")
        g.to_geotiff(fpath, compress=None)
        gnew = karta.read_geotiff(fpath, in_memory=False)

        self.assertEqual(g.transform, gnew.transform)
        self.assertEqual(g.values.dtype, gnew.values.dtype)
        self.assertEqual(g.size, gnew.size)
        self.assertTrue(np.all(g[10:50:3, 15:45:2] == gnew[10:50:3, 15:45:2]))
        self.assertTrue(np.all(g[10:50:3, 45:15:-2] == gnew[10:50:3, 45:15:-2]))
        self.assertTrue(np.all(g[:,:] == gnew[:,:]))
        return

    def test_write_compress(self):
        v = peaks(500)[:100,:]
        utm7 = karta.crs.ProjectedCRS("+proj=utm +zone=7 +north +datum=WGS84",
                                      "UTM 7N (WGS 84)")
        g = karta.RegularGrid([15.0, 15.0, 30.0, 30.0, 0.0, 0.0], v, crs=utm7)

        fpath = os.path.join(TMPDATA, "test.tif")
        g.to_geotiff(fpath, compress="LZW")
        self.assertTrue(os.path.isfile(fpath))
        os.remove(fpath)
        g.to_geotiff(fpath, compress="PACKBITS")
        self.assertTrue(os.path.isfile(fpath))
        return

    def test_read_as_bands(self):
        # write several files and then read as a single multiband grid
        v = peaks(500)[:100,:]
        utm7 = karta.crs.ProjectedCRS("+proj=utm +zone=7 +north +datum=WGS84",
                                      "UTM 7N (WGS 84)")
        g1 = karta.RegularGrid([15.0, 15.0, 30.0, 30.0, 0.0, 0.0], v, crs=utm7)
        g2 = karta.RegularGrid([15.0, 15.0, 30.0, 30.0, 0.0, 0.0], v**2, crs=utm7)
        g3 = karta.RegularGrid([15.0, 15.0, 30.0, 30.0, 0.0, 0.0], v+2, crs=utm7)
        g4 = karta.RegularGrid([15.0, 15.0, 30.0, 30.0, 0.0, 0.0], v*2, crs=utm7)

        paths = []
        for i, g in enumerate((g1, g2, g3, g4)):
            fpath = os.path.join(TMPDATA, "test{0}.tif".format(i))
            g.to_geotiff(fpath, compress=None)
            paths.append(fpath)

        gnew = karta.from_geotiffs(*paths)
        self.assertTrue("+proj=utm" in gnew.crs.get_proj4())
        self.assertTrue("+zone=7" in gnew.crs.get_proj4())
        self.assertEqual(g.transform, gnew.transform)
        self.assertEqual(g.values.dtype, gnew.values.dtype)
        self.assertEqual(gnew.size, (100, 500))
        self.assertEqual(gnew.nbands, 4)
        return


class GdalVirtualArrayTests(unittest.TestCase):

    def setUp(self):
        v = peaks(500)[:100,:]
        utm7 = karta.crs.ProjectedCRS("+proj=utm +zone=7 +north +datum=WGS84",
                                      "UTM 7N (WGS 84)")
        g = karta.RegularGrid([15.0, 15.0, 30.0, 30.0, 0.0, 0.0], v, crs=utm7)

        fpath = os.path.join(TMPDATA, "test.tif")
        g.to_geotiff(fpath, compress=None)
        self.grid = karta.read_geotiff(fpath, in_memory=False)

    def test_slicing_virtual(self):
        """ make sure that slicing a disk-based array works """
        self.grid[5:10, 7:15]
        self.grid[5:10:2, 7:15]
        self.grid[10:5:-1, 7:15]

        a1 = self.grid[12, 7:15, 0]
        self.assertEqual(a1.shape, (8,))
        a2 = self.grid[5:10, 9, 0]
        self.assertEqual(a2.shape, (5,))

        b = self.grid[12, 13, 0]
        self.assertEqual(type(b), np.float64)
        return

    def test_iteration_virtual(self):
        i = 0
        for row in self.grid.values:
            i += 1
        self.assertEqual(i, 100)
        return

def peaks(n=49):
    """ 2d peaks function of MATLAB logo fame. """
    X, Y = np.meshgrid(np.linspace(-3, 3, n), np.linspace(-3, 3, n))
    return 3.0 * (1-X)**2 * np.exp(-X**2 - (Y+1)**2) \
            - 10.0 * (X/5.0 - X**3 - Y**5) * np.exp(-X**2 - Y**2) \
            - 1.0/3.0 * np.exp(-(X+1)**2 - Y**2)

if __name__ == "__main__":
    unittest.main()
