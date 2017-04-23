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

    #def test_setblock_getblock_striped(self):
    #    x, y = np.meshgrid(np.arange(832), np.arange(1024))
    #    d = (x**2+np.sqrt(y))[::2, ::3]

    #    band = self.type((1024, 1024), np.float64, **self.initkwargs)
    #    bi = BandIndexer([band])
    #    bi[::2, 128:960:3, 0] = d

    #    self.assertEqual(np.sum(bi[::2, 128:960:3, 0]-d), 0.0)
    #    return

    def test_setblock_getblock_full(self):
        x, y = np.meshgrid(np.arange(1024), np.arange(1024))
        d = x**2+np.sqrt(y)

        band = self.type((1024, 1024), np.float64, **self.initkwargs)
        band.setblock(0, 0, d)

        self.assertEqual(np.sum(band.getblock(0, 0, 1024, 1024) - d), 0.0)
        return

    def test_setblock_getblock_partial(self):
        x, y = np.meshgrid(np.arange(1024), np.arange(832))
        d = x**2+np.sqrt(y)

        band = self.type((1024, 1024), np.float64, **self.initkwargs)
        band.setblock(128, 0, d)

        self.assertEqual(np.sum(band.getblock(128, 0, 832, 1024) - d), 0.0)
        return

    def test_get_scalar(self):
        x, y = np.meshgrid(np.arange(1024), np.arange(1024))
        d = x**2+np.sqrt(y)

        band = self.type((1024, 1024), np.float64, **self.initkwargs)
        band.setblock(0, 0, d)

        self.assertTrue(band.getblock(4, 3, 1, 1).shape, (1, 1))
        self.assertEqual(band.getblock(4, 3, 1, 1)[0], 11.0)
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

    def test_get_set_typeerror(self):
        bands = [CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32)]
        indexer = BandIndexer(bands)

        with self.assertRaises(TypeError):
            indexer[None]

        with self.assertRaises(TypeError):
            indexer[None] = 3.0

    def test_get_masked(self):
        values = np.ones([16, 16])
        band = CompressedBand((16, 16), np.float32)
        band.setblock(0, 0, values)
        indexer = BandIndexer([band])

        mask = np.zeros([16, 16], dtype=np.bool)
        mask[8:, 2:] = True

        self.assertEqual(np.sum(indexer[mask]), 112)

    def test_set_masked(self):
        values = np.ones([16, 16])
        band = CompressedBand((16, 16), np.float32)
        band.setblock(0, 0, values)
        indexer = BandIndexer([band])

        mask = np.zeros([16, 16], dtype=np.bool)
        mask[8:, 2:] = True

        indexer[mask] = -1
        self.assertEqual(np.sum(indexer[:,:]), 32)

    def test_get_multibanded_errors(self):
        values = np.ones([16, 16])
        bands = [CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32)]
        indexer = BandIndexer(bands)
        with self.assertRaises(TypeError):
            indexer["a", 0, 0]
        with self.assertRaises(TypeError):
            indexer[0, [1, 3], 0]
        with self.assertRaises(TypeError):
            indexer[0, 0, None]

    def test_set_multibanded_errors(self):
        values = np.ones([16, 16])
        bands = [CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32)]
        indexer = BandIndexer(bands)
        with self.assertRaises(TypeError):
            indexer["a", 0, 0] = 1.0
        with self.assertRaises(TypeError):
            indexer[0, [1, 3], 0] = 1.0
        with self.assertRaises(TypeError):
            indexer[0, 0, None] = 1.0

    def test_get_multibanded_integer(self):
        values = np.ones([16, 16])
        bands = [CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32)]
        bands[0].setblock(0, 0, values)
        bands[1].setblock(0, 0, 2*values)
        bands[2].setblock(0, 0, 3*values)
        indexer = BandIndexer(bands)

        result = indexer[10,2]
        self.assertTrue(np.all(result == np.array([1.0, 2.0, 3.0])))

        result = indexer[10,2,:]
        self.assertTrue(np.all(result == np.array([1.0, 2.0, 3.0])))

        # make sure it works with a scalar band index
        result = indexer[10,2,1]
        self.assertEqual(result, 2.0)
        return

    def test_get_multibanded_slice(self):
        values = np.ones([16, 16])
        bands = [CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32)]
        bands[0].setblock(0, 0, values)
        bands[1].setblock(0, 0, 2*values)
        bands[2].setblock(0, 0, 3*values)
        indexer = BandIndexer(bands)

        result = indexer[4:7,2:8,:]
        self.assertEqual(result.shape, (3, 6, 3))
        self.assertTrue(np.all(result[0,0,:] == np.array([1.0, 2.0, 3.0])))

        result = indexer[12:]
        self.assertEqual(result.shape, (4, 16, 3))
        self.assertTrue(np.all(result[0,0,:] == np.array([1.0, 2.0, 3.0])))

        # make sure it works with a scalar band index
        result = indexer[4:7,2:8,1]
        self.assertEqual(result.shape, (3, 6))
        self.assertTrue(np.all(result == 2.0))
        return

    def test_set_multibanded_integer(self):
        values = np.ones([16, 16])
        bands = [CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32)]
        indexer = BandIndexer(bands)

        indexer[3,4,0] = -1.0
        self.assertEqual(bands[0].getblock(3, 4, 1, 1)[0], -1.0)

        indexer[4,5,:] = np.array([1.0, 2.0, 3.0])
        self.assertEqual(bands[0].getblock(4, 5, 1, 1)[0], 1.0)
        self.assertEqual(bands[1].getblock(4, 5, 1, 1)[0], 2.0)
        self.assertEqual(bands[2].getblock(4, 5, 1, 1)[0], 3.0)

        indexer[5,6,:] = 5.0
        self.assertEqual(bands[0].getblock(5, 6, 1, 1)[0], 5.0)
        self.assertEqual(bands[1].getblock(5, 6, 1, 1)[0], 5.0)
        self.assertEqual(bands[2].getblock(5, 6, 1, 1)[0], 5.0)

    def test_set_multibanded_slice(self):
        values = np.ones([16, 16])
        bands = [CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32)]
        indexer = BandIndexer(bands)

        indexer[:,:,0] = values
        indexer[:,:,1:] = 2*values
        self.assertTrue(np.all(bands[0].getblock(0, 0, 16, 16) == 1.0))
        self.assertTrue(np.all(bands[1].getblock(0, 0, 16, 16) == 2.0))
        self.assertTrue(np.all(bands[2].getblock(0, 0, 16, 16) == 2.0))

        indexer[:,:,1:] = np.dstack([2*values, 3*values])
        self.assertTrue(np.all(bands[1].getblock(0, 0, 16, 16) == 2.0))
        self.assertTrue(np.all(bands[2].getblock(0, 0, 16, 16) == 3.0))
        return

    def test_set_multibanded_masked_scalar(self):
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

    def test_set_multibanded_masked_array(self):
        values = np.ones([16, 16])
        bands = [CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32)]

        mask = np.zeros([16, 16], dtype=np.bool)
        mask[8:, 2:] = True

        indexer = BandIndexer(bands)
        indexer[:,:,:] = np.zeros([16, 16])
        indexer[mask] = np.ones(8*14)

        self.assertEqual(np.sum(indexer[:,:]), 336)
        return

    def test_shape(self):
        bands = [CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32),
                 CompressedBand((16, 16), np.float32)]

        indexer1 = BandIndexer([bands[0]])
        self.assertEqual(indexer1.shape, (16, 16))

        indexer3 = BandIndexer(bands)
        self.assertEqual(indexer3.shape, (16, 16))
        return

if __name__ == "__main__":
    unittest.main()
