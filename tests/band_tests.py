import unittest
import numpy as np

from karta.raster import SimpleBand, CompressedBand
from karta.raster.band import BandIndexer

class GenericBandTests(object):
    """ Tests that all Band classes must pass """

    def test_get_dtype(self):
        band = self.type((64, 64), np.float64, **self.initkwargs)
        self.assertEqual(band.dtype, np.float64)
        return

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

class BandIndexerTests(unittest.TestCase):

    def test_get_masked(self):
        values = np.ones([16, 16])
        band = CompressedBand((16, 16), np.float32)
        band[:,:] = values
        indexer = BandIndexer([band])

        mask = np.zeros([16, 16], dtype=np.bool)
        mask[8:, 2:] = True

        self.assertEqual(np.sum(indexer[mask]), 112)

    def test_set_masked(self):
        values = np.ones([16, 16])
        band = CompressedBand((16, 16), np.float32)
        band[:,:] = values
        indexer = BandIndexer([band])

        mask = np.zeros([16, 16], dtype=np.bool)
        mask[8:, 2:] = True

        indexer[mask] = -1
        self.assertEqual(np.sum(indexer[:,:]), 32)

if __name__ == "__main__":
    unittest.main()
