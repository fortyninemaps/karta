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

    def test_get_scalar(self):

        x, y = np.meshgrid(np.arange(1024), np.arange(1024))
        d = x**2+np.sqrt(y)

        band = self.type((1024, 1024), np.float64, **self.initkwargs)
        band[:, :] = d

        self.assertEqual(band[4,3], 11.0)
        self.assertTrue(band[4:5,3:4].shape, (1, 1))
        self.assertEqual(band[-1,-1], d[-1,-1])
        return

    def test_initval(self):
        band = self.type((1024, 1024), np.float64, initval=0.0)
        self.assertTrue(band is not None)
        return


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

    def test_get_multibanded(self):
        values = np.ones([16, 16])
        bands = [CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32)]
        bands[0][:,:] = values
        bands[1][:,:] = 2*values
        bands[2][:,:] = 3*values

        indexer = BandIndexer(bands)
        result = indexer[4:7,2:8,:]
        self.assertEqual(result.shape, (3, 6, 3))
        self.assertTrue(np.all(result[0,0,:] == np.array([1.0, 2.0, 3.0])))

        # make sure it works with a scalar band index
        result = indexer[4:7,2:8,1]
        self.assertEqual(result.shape, (3, 6))
        self.assertTrue(np.all(result == 2.0))
        return

    def test_set_multibanded(self):
        values = np.ones([16, 16])
        bands = [CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32)]

        indexer = BandIndexer(bands)
        indexer[:,:,0] = values
        indexer[:,:,1:] = 2*values

        self.assertTrue(np.all(bands[0][:,:] == 1.0))
        self.assertTrue(np.all(bands[1][:,:] == 2.0))
        self.assertTrue(np.all(bands[2][:,:] == 2.0))

        indexer[:,:,1:] = np.dstack([2*values, 3*values])

        self.assertTrue(np.all(bands[1][:,:] == 2.0))
        self.assertTrue(np.all(bands[2][:,:] == 3.0))
        return

    def test_set_multibanded_masked(self):
        values = np.ones([16, 16])
        bands = [CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32)]

        mask = np.zeros([16, 16], dtype=np.bool)
        mask[8:, 2:] = True

        indexer = BandIndexer(bands)
        indexer[:,:] = np.zeros([16, 16])
        indexer[mask] = 1.0

        self.assertEqual(np.sum(indexer[:,:]), 336)
        return

if __name__ == "__main__":
    unittest.main()
