import unittest
import os.path
import numpy as np
from test_helper import TESTDATA

import karta
from karta.raster import _gtiff

class GdalTests(unittest.TestCase):

    def test_numpy_type_coercion(self):
        self.assertEqual(_gtiff.numpy_dtype(2), np.uint16)
        self.assertEqual(_gtiff.numpy_dtype(3), np.int16)
        self.assertEqual(_gtiff.numpy_dtype(4), np.uint32)
        self.assertEqual(_gtiff.numpy_dtype(5), np.int32)
        self.assertEqual(_gtiff.numpy_dtype(6), np.float32)
        self.assertEqual(_gtiff.numpy_dtype(7), np.float64)
        self.assertEqual(_gtiff.numpy_dtype(8), np.complex64)
        self.assertEqual(_gtiff.numpy_dtype(9), np.complex64)
        self.assertEqual(_gtiff.numpy_dtype(10), np.complex64)
        self.assertEqual(_gtiff.numpy_dtype(11), np.complex64)
        return

    def test_io(self):
        # try writing a file, then read it back in and verify that it matches
        v = karta.raster.peaks(500)[:100,:]
        utm7 = karta.crs.Proj4CRS("+proj=utm +zone=7 +north", "+ellps=WGS84")
        g = karta.RegularGrid([15.0, 15.0, 30.0, 30.0, 0.0, 0.0], v, crs=utm7)

        fpath = os.path.join(TESTDATA, "test.tif")
        g.to_gtiff(fpath, compress=None)
        gnew = karta.read_gtiff(fpath)

        self.assertTrue("+proj=utm" in gnew.crs.project.srs)
        self.assertTrue("+zone=7" in gnew.crs.project.srs)
        self.assertEqual(g.transform, gnew.transform)
        self.assertEqual(g.values.dtype, gnew.values.dtype)
        self.assertTrue(np.all(g.values == gnew.values))
        return

    def test_io_virtual(self):
        # try writing a file, then open it without loading into memory and verify
        v = karta.raster.peaks(500)[:100,:]
        utm7 = karta.crs.Proj4CRS("+proj=utm +zone=7 +north", "+ellps=WGS84")
        g = karta.RegularGrid([15.0, 15.0, 30.0, 30.0, 0.0, 0.0], v, crs=utm7)

        fpath = os.path.join(TESTDATA, "test.tif")
        g.to_gtiff(fpath, compress=None)
        gnew = karta.read_gtiff(fpath, in_memory=False)

        self.assertEqual(g.transform, gnew.transform)
        self.assertEqual(g.values.dtype, gnew.values.dtype)
        self.assertEqual(g.size, gnew.size)
        self.assertTrue(np.all(g.values[10:50:3, 15:45:2] ==
                               gnew.values[10:50:3, 15:45:2]))
        self.assertTrue(np.all(g.values[10:50:3, 45:15:-2] ==
                               gnew.values[10:50:3, 45:15:-2]))
        self.assertTrue(np.all(g.values[:,:] == gnew.values[:,:]))
        return

    def test_write_compress(self):
        v = karta.raster.peaks(500)[:100,:]
        utm7 = karta.crs.Proj4CRS("+proj=utm +zone=7 +north", "+ellps=WGS84")
        g = karta.RegularGrid([15.0, 15.0, 30.0, 30.0, 0.0, 0.0], v, crs=utm7)

        fpath = os.path.join(TESTDATA, "test.tif")
        g.to_gtiff(fpath, compress="LZW")
        g.to_gtiff(fpath, compress="PACKBITS")
        return

if __name__ == "__main__":
    unittest.main()
