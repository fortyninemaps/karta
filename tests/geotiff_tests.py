import unittest
import numpy as np
from test_helper import TESTDATA

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

