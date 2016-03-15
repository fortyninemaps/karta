import unittest
import numpy as np

from karta.raster import SimpleBand, CompressedBand

class GenericBandTests(object):
    """ Tests that all Band classes must pass """

    def test_setblock_getblock_full(self):

        x, y = np.meshgrid(np.arange(1024), np.arange(1024))
        d = x**2+np.sqrt(y)

        band = self.type((1024, 1024), np.float64, **self.initkwargs)
        band[:, :] = d

        self.assertEqual(np.sum(band[:,:] - d), 0.0)

    def test_setblock_getblock_partial(self):

        x, y = np.meshgrid(np.arange(1024), np.arange(832))
        d = x**2+np.sqrt(y)

        band = self.type((1024, 1024), np.float64, **self.initkwargs)
        band[128:960, :] = d

        self.assertEqual(np.sum(band[128:960,:]-d), 0.0)

    def test_setblock_getblock_striped(self):

        x, y = np.meshgrid(np.arange(832), np.arange(1024))
        d = (x**2+np.sqrt(y))[::2, ::3]

        band = self.type((1024, 1024), np.float64, **self.initkwargs)
        band[::2, 128:960:3] = d

        self.assertEqual(np.sum(band[::2,128:960:3]-d), 0.0)



class SimpleBandTests(unittest.TestCase, GenericBandTests):

    def setUp(self):
        self.type = SimpleBand
        self.initkwargs = dict()


class CompressedBandTests(unittest.TestCase, GenericBandTests):

    def setUp(self):
        self.type = CompressedBand
        self.initkwargs = dict(chunksize=(256, 256))

if __name__ == "__main__":
    unittest.main()
