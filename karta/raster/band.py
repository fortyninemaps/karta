import blosc
import numpy as np
from math import ceil

class SimpleBand(object):
    """ SimpleBand wraps a numpy.ndarray for storage. """

    def __init__(self, size, dtype):
        self.size = size
        self.array = np.empty(size, dtype=dtype)
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

    def __init__(self, size, dtype, chunksize=(256, 256)):
        assert len(size) == 2
        self.size = size
        self.dtype = dtype
        self._chunksize = chunksize

        self.nchunkrows = int(ceil(float(size[0])/float(chunksize[0])))
        self.nchunkcols = int(ceil(float(size[1])/float(chunksize[1])))
        nchunks = self.nchunkrows * self.nchunkcols

        # Data store
        self._data = [None for i in range(nchunks)]

        # 0 => unset
        # 1 => set
        self.chunkstatus = np.zeros(nchunks, dtype=np.int8)

    def __getitem__(self, key):

        if isinstance(key, int):
            irow = key // self.size[1]
            icol = key % self.size[1]
            return self._getblock(irow, icol, (1, 1))[0]

        elif isinstance(key, tuple):

            if len(key) != 2:
                raise IndexError("band can only be indexed along two dimensions")

            k0, k1 = key

            if isinstance(k0, int):
                yoff = k0
                ny = 1
                sy = 1

            elif isinstance(k0, slice):
                if k0.start is None:
                    yoff = 0
                else:
                    yoff = k0.start
                if k0.stop is None:
                    ny = self.size[0]-yoff
                else:
                    ny = k0.stop-yoff
                if k0.step is None:
                    sy = 1
                else:
                    sy = k0.step

            else:
                raise IndexError("slicing with instances of '{0}' not "
                                 "supported".format(type(k0)))

            if isinstance(k1, int):
                xoff = k1
                nx = 1
                sx = 1

            elif isinstance(k1, slice):
                if k1.start is None:
                    xoff = 0
                else:
                    xoff = k1.start
                if k1.stop is None:
                    nx = self.size[1]-xoff
                else:
                    nx = k1.stop-xoff
                if k1.step is None:
                    sx = 1
                else:
                    sx = k1.step

            else:
                raise IndexError("slicing with instances of '{0}' not "
                                 "supported".format(type(k1)))

            return self._getblock(yoff, xoff, (ny, nx))[::sy,::sx]

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

            k0, k1 = key

            if isinstance(k0, int):
                yoff = k0
                ny = 1
                sy = 1

            elif isinstance(k0, slice):
                if k0.start is None:
                    yoff = 0
                else:
                    yoff = k0.start
                if k0.stop is None:
                    ny = self.size[0]-yoff
                else:
                    ny = k0.stop-yoff
                if k0.step is None:
                    sy = 1
                else:
                    sy = k0.step

            else:
                raise IndexError("slicing with instances of '{0}' not "
                                 "supported".format(type(k0)))

            if isinstance(k1, int):
                xoff = k1
                nx = 1
                sx = 1

            elif isinstance(k1, slice):
                if k1.start is None:
                    xoff = 0
                else:
                    xoff = k1.start
                if k1.stop is None:
                    nx = self.size[1]-xoff
                else:
                    nx = k1.stop-xoff
                if k1.step is None:
                    sx = 1
                else:
                    sx = k1.step

            else:
                raise IndexError("slicing with instances of '{0}' not "
                                 "supported".format(type(k1)))

            # TODO add size checks

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
        size = array.shape
        chunksize = self._chunksize
        chunkrowstart = yoff // chunksize[0]
        chunkcolstart = xoff // chunksize[1]

        for i, yst, yen, xst, xen in self._getchunks(yoff, xoff, *size):

            # Get from data store
            if self.chunkstatus[i] != self.CHUNKUNSET:
                chunkdata = self._retrieve(i)
            else:
                chunkdata = np.zeros(self._chunksize, dtype=self.dtype)

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
        result = np.empty(size, self.dtype)
        for i, yst, yen, xst, xen in self._getchunks(yoff, xoff, *size):

            # Compute the bounds in the output
            oy0 = max(0, yst-yoff)
            oy1 = min(size[0], yen-yoff)
            ox0 = max(0, xst-xoff)
            ox1 = min(size[1], xen-xoff)

            if self.chunkstatus[i] == self.CHUNKUNSET:
                result[oy0:oy1, ox0:ox1] = np.zeros((oy1-oy0, ox1-ox0),
                                                    dtype=self.dtype)

            else:
                # Compute the extents from the chunk to retain
                cy0 = max(yoff, yst) - yst
                cy1 = min(yoff+size[0], yen) - yst
                cx0 = max(xoff, xst) - xst
                cx1 = min(xoff+size[1], xen) - xst

                result[oy0:oy1, ox0:ox1] = self._retrieve(i)[cy0:cy1, cx0:cx1]

        return result

