"""
Band implementations for storing data in Karta Grid instances

Overview
--------

`BandIndexer` interface for accessing data from one or more bands

`SimpleBand` use numpy arrays for data storage

`CompressedBand` uses blosc compression to reduce in-memory footprint

Implementation
--------------

Bands are expected to implement the following methods:
    - `__init__(self, size, dtype, initval=None)`
    - `getblock(self, yoff, xoff, ny, nx)`
    - `setblock(self, yoff, xoff, array)`

Attributes:
    - `dtype`
    - `size`

The following methods are deprecated:
    - `__getitem__(self, key)`, accepting as *key* any of
        - an int
        - a slice
        - a 2-tuple of ints
        - a 2-tuple of slices
    - `__setitem__(self, key, value)`, accepting as *key* the same
      possibilities as __getitem__
"""

import blosc
import numpy as np
from numbers import Real, Integral
from math import ceil

class BandIndexer(object):

    def __init__(self, bands):
        self.bands = bands

    def __getitem__(self, key):
        if isinstance(key, np.ndarray):
            return self._get_from_array_mask(key)

        if isinstance(key, slice):
            key = (key, slice(None, None, None), slice(None, None, None))

        if not isinstance(key, tuple):
            raise TypeError("key should be an array or a tuple")

        collapse_rows = collapse_cols = collapse_bands = False
        ny, nx = self.bands[0].size

        if isinstance(key[0], Integral):
            collapse_rows = True
            r = key[0] % ny
            ystart, yend, ystep = (r, r+1, 1)
        elif isinstance(key[0], slice):
            ystart, yend, ystep = key[0].indices(ny)
        else:
            raise TypeError("first key item should be an integer or a slice")

        if isinstance(key[1], Integral):
            collapse_cols = True
            r = key[1] % nx
            xstart, xend, xstep = (r, r+1, 1)
        elif isinstance(key[1], slice):
            xstart, xend, xstep = key[1].indices(nx)
        else:
            raise TypeError("second key item should be an integer or a slice")

        if len(key) == 2:
            bands = list(range(len(self.bands)))
        elif len(key) == 3 and isinstance(key[2], Integral):
            collapse_bands = True
            bands = [key[2] % len(self.bands)]
        elif len(key) == 3 and isinstance(key[2], slice):
            bands = list(range(*key[2].indices(len(self.bands))))
        else:
            raise TypeError("third key item should be an integer or a slice")

        if ystep < 0:
            ystart, yend = yend+1, ystart+1

        if xstep < 0:
            xstart, xend = xend+1, xstart+1

        shape = [1 + (yend-ystart-1) // abs(ystep),
                 1 + (xend-xstart-1) // abs(xstep),
                 len(bands)]

        out = np.empty(shape, dtype = self.bands[0].dtype)

        for i, iband in enumerate(bands):
            band = self.bands[iband]
            band_values = band.getblock(ystart, xstart, yend-ystart, xend-xstart)
            out[:,:,i] = band_values[::ystep,::xstep]

        if collapse_bands:
            out = out[:,:,0]
        if collapse_cols:
            out = out[:,0]
        if collapse_rows:
            out = out[0]

        return out

    def __setitem__(self, key, value):
        if isinstance(key, np.ndarray):
            return self._set_from_array_mask(key, value)

        if isinstance(key, slice):
            key = (key, slice(None, None, None), slice(None, None, None))

        if not isinstance(key, tuple):
            raise TypeError("key should be an array or a tuple")

        ny, nx = self.bands[0].size

        if isinstance(key[0], Integral):
            r = key[0] % ny
            ystart, yend, ystep = (r, r+1, 1)
        elif isinstance(key[0], slice):
            ystart, yend, ystep = key[0].indices(ny)
        else:
            raise TypeError("first key item should be an integer or a slice")

        if isinstance(key[1], Integral):
            r = key[1] % nx
            xstart, xend, xstep = (r, r+1, 1)
        elif isinstance(key[1], slice):
            xstart, xend, xstep = key[1].indices(nx)
        else:
            raise TypeError("second key item should be an integer or a slice")

        if len(key) == 2:
            bands = list(range(len(self.bands)))
        elif len(key) == 3 and isinstance(key[2], Integral):
            collapse_bands = True
            bands = [key[2] % len(self.bands)]
        elif len(key) == 3 and isinstance(key[2], slice):
            bands = list(range(*key[2].indices(len(self.bands))))
        else:
            raise TypeError("third key item should be an integer or a slice")

        if not (xstep == ystep == 1):
            raise NotImplementedError("setting band values with stepped slices")

            #if ystep < 0:
            #    ystart, yend = yend+1, ystart+1

            #if xstep < 0:
            #    xstart, xend = xend+1, xstart+1

        shape = [1 + (yend-ystart-1) // abs(ystep),
                 1 + (xend-xstart-1) // abs(xstep),
                 len(bands)]

        if isinstance(value, np.ndarray) and (value.ndim == 1) and (shape[0] == shape[1] == 1):
            val_array = np.reshape(np.atleast_3d(value), shape)
        else:
            val_array = np.broadcast_to(np.atleast_3d(value), shape)

        for i, iband in enumerate(bands):
            band = self.bands[iband]
            band.setblock(ystart, xstart, val_array[:,:,i])

        return

    def _get_from_array_mask(self, mask):
        # The mask is assumed to be in (row, column[, band]) order
        # TODO: make this memory efficient
        if mask.ndim == 2:
            return self[:,:,:][mask]
        elif mask.ndim == 3:
            return self[:,:,:][mask]
        else:
            raise IndexError("masking array must have two or three dimensions")

    def _set_from_array_mask(self, mask, value):
        # The mask is assumed to be in (row, column[, band]) order
        # TODO: make this memory efficient
        for i, band in enumerate(self.bands):

            if mask.ndim == 3:
                mask_ = mask[:,:,i]
            else:
                mask_ = mask

            tmp = band.getblock(0, 0, *band.size)
            if isinstance(value, Real) or (value.ndim == 1):
                tmp[mask_] = value
            else:
                tmp[mask_] = value[:,i]
            band.setblock(0, 0, tmp)

    def __iter__(self):
        nx = self.bands[0].size[1]
        for i in range(self.bands[0].size[0]):
            if len(self.bands) == 1:
                yield self.bands[0].getblock(i, 0, 1, nx)
            else:
                yield np.vstack([b.getblock(i, 0, 1, nx) for b in self.bands])

    @property
    def shape(self):
        """ Returns the dimensions of raster bands. If there is a single
        (m x n) band, output is a tuple (m, n). If there are N>1 bands, output
        is a tuple (N, m, n).
        """
        if len(self.bands) == 0:
            raise ValueError("no bands")
        else:
            return self.bands[0].size

    @property
    def dtype(self):
        """ Returns bands' dtype """
        return self.bands[0].dtype

class SimpleBand(object):
    """ SimpleBand wraps a numpy.ndarray for storage. """

    def __init__(self, size, dtype, initval=None):
        self.size = size
        if initval is None:
            self._array = np.empty(size, dtype=dtype)
        else:
            self._array = np.full(size, initval, dtype=dtype)
        self.dtype = dtype

    def getblock(self, yoff, xoff, ny, nx):
        return self._array[yoff:yoff+ny, xoff:xoff+nx]

    def setblock(self, yoff, xoff, array):
        (ny, nx) = array.shape
        self._array[yoff:yoff+ny, xoff:xoff+nx] = array
        return

class CompressedBand(object):
    """ CompressedBand is a chunked, blosc-compressed array. """
    CHUNKSET = 1
    CHUNKUNSET = 0

    def __init__(self, size, dtype, chunksize=(256, 256), initval=0):
        """ Initialize a CompressedBand instance.

        Parameters
        ----------
        size : tuple of two ints
            size of band in pixels
        dtype : type
            data type of pixel values
        chunksize : tuple of two ints, optional
            size of compressed chunks, default (256, 256)
        initval : value, optional
            if set, the entire grid is initialized with this value, which should
            be of *dtype*
        """
        assert len(size) == 2
        self.size = size
        self.dtype = dtype
        self._chunksize = chunksize
        self._initval = initval

        self.nchunkrows = int(ceil(float(size[0])/float(chunksize[0])))
        self.nchunkcols = int(ceil(float(size[1])/float(chunksize[1])))
        nchunks = self.nchunkrows * self.nchunkcols

        # Data store
        self._data = [None for i in range(nchunks)]

        # 0 => unset
        # 1 => set
        self.chunkstatus = np.zeros(nchunks, dtype=np.int8)
        return

    def _store(self, array, index):
        self._data[index] = blosc.compress(array.tostring(),
                                           np.dtype(self.dtype).itemsize)
        self.chunkstatus[index] = self.CHUNKSET
        return

    def _retrieve(self, index):
        bytestr = blosc.decompress(self._data[index])
        return np.fromstring(bytestr, dtype=self.dtype).reshape(self._chunksize)

    def _getchunks(self, yoff, xoff, ny, nx):
        """ Return a generator returning tuples identifying chunks covered by a
        range. The tuples contain (chunk_number, ystart, yend, xstart, xend)
        for each chunk touched by a region defined by corner indices and region
        size. """
        chunksize = self._chunksize
        ystart = yoff // chunksize[0]
        yend = ceil(float(yoff+ny) / chunksize[0])
        xstart = xoff // chunksize[1]
        xend = ceil(float(xoff+nx) / chunksize[1])

        nxchunks = int(ceil(float(self.size[1])/float(chunksize[1])))

        i = ystart
        while i < yend:

            j = xstart

            while j < xend:
                chunk_number = i*nxchunks + j
                chunk_ystart = i*chunksize[0]
                chunk_xstart = j*chunksize[1]
                chunk_yend = min((i+1)*chunksize[0], self.size[0])
                chunk_xend = min((j+1)*chunksize[1], self.size[1])
                yield (chunk_number, chunk_ystart, chunk_yend,
                       chunk_xstart, chunk_xend)
                j += 1

            i+= 1

    def setblock(self, yoff, xoff, array):
        """ Store block of values in *array* starting at offset *yoff*, *xoff*.
        """
        size = array.shape[:2]
        chunksize = self._chunksize

        for i, yst, yen, xst, xen in self._getchunks(yoff, xoff, *size):

            # Get from data store
            if self.chunkstatus[i] == self.CHUNKSET:
                chunkdata = self._retrieve(i)
            else:
                chunkdata = np.full(self._chunksize, self._initval, dtype=self.dtype)

            # Compute region within chunk to place data in
            cy0 = max(0, yoff-yst)
            cy1 = min(chunksize[0], yoff+size[0]-yst)
            cx0 = max(0, xoff-xst)
            cx1 = min(chunksize[1], xoff+size[1]-xst)

            # Compute region to slice from data
            dy0 = max(0, yst-yoff)
            dy1 = min(size[0], yen-yoff)
            dx0 = max(0, xst-xoff)
            dx1 = min(size[1], xen-xoff)

            chunkdata[cy0:cy1, cx0:cx1] = array[dy0:dy1, dx0:dx1]

            # Return to data store
            self._store(chunkdata, i)
        return

    def getblock(self, yoff, xoff, ny, nx):
        """ Retrieve values with dimensions *size*, starting at offset *yoff*,
        *xoff*.
        """
        result = np.empty([ny, nx], self.dtype)
        for i, yst, yen, xst, xen in self._getchunks(yoff, xoff, ny, nx):

            # Compute the bounds in the output
            oy0 = max(0, yst-yoff)
            oy1 = min(ny, yen-yoff)
            ox0 = max(0, xst-xoff)
            ox1 = min(nx, xen-xoff)

            if self.chunkstatus[i] == self.CHUNKUNSET:
                result[oy0:oy1, ox0:ox1] = np.full((oy1-oy0, ox1-ox0),
                                                   self._initval,
                                                   dtype=self.dtype)

            else:
                # Compute the extents from the chunk to retain
                cy0 = max(yoff, yst) - yst
                cy1 = min(yoff+ny, yen) - yst
                cx0 = max(xoff, xst) - xst
                cx1 = min(xoff+nx, xen) - xst

                result[oy0:oy1, ox0:ox1] = self._retrieve(i)[cy0:cy1, cx0:cx1]

        return result

