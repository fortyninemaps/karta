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
    - `__getitem__(self, key)`, accepting as *key* any of
        - an int
        - a slice
        - a 2-tuple of ints
        - a 2-tuple of slices
    - `__setitem__(self, key, value)`, accepting as *key* the same
      possibilities as __getitem__

Additionally, a Band should have a `dtype` attribute indicating the type of the
stored data.
"""

import blosc
import numpy as np
from numbers import Integral
from math import ceil

class BandIndexer(object):

    def __init__(self, bands):
        self.bands = bands

    def __getitem__(self, key):
        if len(self.bands) == 1:
            if isinstance(key, np.ndarray):
                return self.bands[0][:,:][key]
            elif isinstance(key, tuple):
                if (len(key) == 3) and (key[2] != 0):
                    raise IndexError("third element of key exceeds number of bands")
                return self.bands[0][key[:2]]
            else:
                return self.bands[0][key]
        else:
            if isinstance(key, np.ndarray):
                return np.dstack([b[:,:][key] for b in self.bands])
            else:
                if len(key) not in (2, 3):
                    raise IndexError("indexing tuple must have length 2 or 3")
                sr, sc = key[:2]
                sb = slice(None, None, None) if len(key) == 2 else key[2]
                if isinstance(sb, slice):
                    return np.dstack([b[sr,sc] for b in self.bands[sb]])
                else:
                    return self.bands[sb][sr,sc]

    def __setitem__(self, key, value):
        if len(self.bands) == 1:
            if isinstance(key, np.ndarray):
                tmp = self.bands[0][:,:]
                tmp[key] = value
                self.bands[0][:,:] = tmp
            else:
                self.bands[0][key] = value
        else:
            if isinstance(key, np.ndarray):
                if isinstance(value, np.ndarray):
                    for b, v in zip(self.bands, value):
                        tmp = b[:,:]
                        tmp[key] = v
                        b[:,:] = tmp
                else:
                    for b in self.bands:
                        tmp = b[:,:]
                        tmp[key] = value
                        b[:,:] = tmp
            else:
                if len(key) not in (2, 3):
                    raise IndexError("indexing tuple must have length 2 or 3")
                sr, sc = key[:2]
                sb = slice(None, None, None) if len(key) == 2 else key[2]
                if isinstance(sb, slice):
                    if len(value.shape) == 3:
                        for ival, iband in enumerate(range(*sb.indices(len(self.bands)))):
                            self.bands[iband][sr,sc] = value[:,:,ival]
                    else:
                        for b in self.bands[sb]:
                            b[sr,sc] = value
                else:
                    self.bands[sb][sr,sc] = value
        return

    def __iter__(self):
        for i in range(self.bands[0].size[0]):
            if len(self.bands) == 1:
                yield self.bands[0][i,:]
            else:
                yield np.vstack([b[i,:] for b in self.bands])

    @property
    def shape(self):
        """ Returns the dimensions of raster bands. If there is a single
        (m x n) band, output is a tuple (m, n). If there are N>1 bands, output
        is a tuple (N, m, n).
        """
        if len(self.bands) == 0:
            raise ValueError("no bands")
        elif len(self.bands) == 1:
            return self.bands[0].size
        else:
            return (len(self.bands), self.bands[0].size[0], self.bands[0].size[1])

    @property
    def dtype(self):
        """ Returns bands' dtype """
        return self.bands[0].dtype

class SimpleBand(object):
    """ SimpleBand wraps a numpy.ndarray for storage. """

    def __init__(self, size, dtype, initval=None):
        self.size = size
        if initval is None:
            self.array = np.empty(size, dtype=dtype)
        else:
            self.array = np.full(size, initval, dtype=dtype)
        self.dtype = dtype

    def __getitem__(self, key):
        return self.array[key]

    def __setitem__(self, key, value):
        self.array[key] = value
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

    def __getitem__(self, key):

        if isinstance(key, Integral):
            irow = key // self.size[1]
            icol = key % self.size[1]
            return self._getblock(irow, icol, (1, 1))[0]

        elif isinstance(key, tuple):

            if len(key) != 2:
                raise IndexError("band can only be indexed along two dimensions")

            kr, kc = key

            if isinstance(kr, Integral):
                yoff = kr % self.size[0]
                ny = 1
                ystride = 1

            elif isinstance(kr, slice):
                ystart, ystop, ystride = kr.indices(self.size[0])
                yoff = min(ystart, ystop)
                ny = abs(ystop-ystart)
                if ystride < 0:
                    yoff += 1

            else:
                raise IndexError("slicing with instances of '{0}' not "
                                 "supported".format(type(kr)))

            if isinstance(kc, Integral):
                xoff = kc % self.size[1]
                nx = 1
                xstride = 1

            elif isinstance(kc, slice):
                xstart, xstop, xstride = kc.indices(self.size[1])
                xoff = min(xstart, xstop)
                nx = abs(xstop-xstart)
                if xstride < 0:
                    xoff += 1

            else:
                raise IndexError("slicing with instances of '{0}' not "
                                 "supported".format(type(kc)))

            if isinstance(kr, Integral) and isinstance(kc, Integral):
                return self._getblock(yoff, xoff, (ny, nx))[0,0]
            else:
                return self._getblock(yoff, xoff, (ny, nx))[::ystride,::xstride]

        elif isinstance(key, slice):
            start, stop, stride = key.indices(self.size[0])
            yoff = min(start, stop)
            if stride < 0:
                yoff += 1
            ny = abs(stop-start)
            return self._getblock(yoff, 0, (ny, self.size[1]))[::stride]

        else:
            raise IndexError("indexing with instances of '{0}' not "
                             "supported".format(type(key)))

    def __setitem__(self, key, value):

        if isinstance(key, int):
            irow = key // self.size[1]
            icol = key % self.size[1]
            self._setblock(irow, icol, np.array(value, dtype=self.dtype))

        elif isinstance(key, tuple):

            if len(key) != 2:
                raise IndexError("band can only be indexed along two dimensions")

            kr, kc = key

            if isinstance(kr, int):
                yoff = kr
                ny = 1
                sy = 1

            elif isinstance(kr, slice):
                yoff, _ystop, sy = kr.indices(self.size[0])
                ny = _ystop - yoff

            else:
                raise IndexError("slicing with instances of '{0}' not "
                                 "supported".format(type(kr)))

            if isinstance(kc, int):
                xoff = kc
                nx = 1
                sx = 1

            elif isinstance(kc, slice):
                xoff, _xstop, sx = kc.indices(self.size[1])
                nx = _xstop - xoff

            else:
                raise IndexError("slicing with instances of '{0}' not "
                                 "supported".format(type(kc)))

            if not hasattr(value, "shape"):
                value = value*np.ones((ny, nx))

            vny, vnx = value.shape[:2]
            if (ceil(float(ny)/sy) != vny) or (ceil(float(nx)/sx) != vnx):
                raise IndexError("Cannot insert array with size ({vny}, {vnx}))"
                        " into slice with size ({ny}, {nx})".format(
                            vny=vny, vnx=vnx, ny=ny, nx=nx))

            if sy == sx == 1:
                self._setblock(yoff, xoff, value)

            else:
                temp = self._getblock(yoff, xoff, (ny, nx))
                temp[::sy,::sx] = value
                self._setblock(yoff, xoff, temp)

        else:
            raise IndexError("indexing with instances of '{0}' not "
                             "supported".format(type(key)))

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

    def _setblock(self, yoff, xoff, array):
        """ Store block of values in *array* starting at offset *yoff*, *xoff*.
        """
        size = array.shape
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

    def _getblock(self, yoff, xoff, size):
        """ Retrieve values with dimensions *size*, starting at offset *yoff*,
        *xoff*.
        """
        result = np.empty(size, self.dtype)
        for i, yst, yen, xst, xen in self._getchunks(yoff, xoff, *size):

            # Compute the bounds in the output
            oy0 = max(0, yst-yoff)
            oy1 = min(size[0], yen-yoff)
            ox0 = max(0, xst-xoff)
            ox1 = min(size[1], xen-xoff)

            if self.chunkstatus[i] == self.CHUNKUNSET:
                result[oy0:oy1, ox0:ox1] = np.full((oy1-oy0, ox1-ox0),
                                                   self._initval,
                                                   dtype=self.dtype)

            else:
                # Compute the extents from the chunk to retain
                cy0 = max(yoff, yst) - yst
                cy1 = min(yoff+size[0], yen) - yst
                cx0 = max(xoff, xst) - xst
                cx1 = min(xoff+size[1], xen) - xst

                result[oy0:oy1, ox0:ox1] = self._retrieve(i)[cy0:cy1, cx0:cx1]

        return result

