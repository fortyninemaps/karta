import blosc
import numpy as np
from math import ceil

class SimpleBand(object):
    """ SimpleBand wraps a numpy.ndarray for storage. """

    def __init__(self, size, dtype):
        self.array = np.empty(size, dtype=dtype)

    def __getitem__(self, key):
        return self.array[key]

    def __setitem__(self, key):
        return self.array[key]

class CompressedBand(object):
    """ CompressedBand is a chunked, blosc-compressed array. """
    CHUNKSET = 1
    CHUNKUNSET = 0

    def __init__(self, size, dtype, chunksize=(256, 256)):
        assert len(size) == 2
        self.size = size
        self.dtype = dtype
        self.chunksize = chunksize

        self.nchunkrows = ceil(float(size[0])/float(chunksize[0]))
        self.nchunkcols = ceil(float(size[1])/float(chunksize[1]))
        nchunks = self.nchunkrows * self.nchunkcols

        # Data store
        self._data = [None for i in range(nchunks)]

        # 0 => unset
        # 1 => set
        self.chunkstatus = np.zeros(nchunks, dtype=np.int8)

    def _store(self, array, index):
        self._data[index] = blosc.compress(array.tostring(), np.dtype(self.dtype).itemsize)
        self.chunkstatus[index] = self.CHUNKSET
        return

    def _retrieve(self, index):
        bytestr = blosc.decompress(self._data[index])
        return np.fromstring(bytestr, dtype=self.dtype).reshape(self.chunksize)

    def getchunks(self, yoff, xoff, ny, nx):
        """ Return a generator returning tuples identifying chunks covered by a
        range. The tuples contain (chunk_number, ystart, yend, xstart, xend)
        for each chunk touched by a region defined by corner indices and region
        size. """
        chunksize = self.chunksize
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
                yield (chunk_number, chunk_ystart, chunk_yend, chunk_xstart, chunk_xend)
                j += 1

            i+= 1

    def setdata(self, yoff, xoff, array):
        size = array.shape
        chunksize = self.chunksize
        chunkrowstart = yoff // chunksize[0]
        chunkcolstart = xoff // chunksize[1]

        for i, yst, yen, xst, xen in self.getchunks(yoff, xoff, *size):

            # Get from data store
            if self.chunkstatus[i] != self.CHUNKUNSET:
                chunkdata = self._retrieve(i)
            else:
                chunkdata = np.zeros(self.chunksize, dtype=self.dtype)

            # Identify chunk
            irow = (i // self.nchunkcols - chunkrowstart)
            icol = (i % self.nchunkcols - chunkcolstart)
            chunk_y0 = irow*chunksize[0]
            chunk_y1 = (irow+1)*chunksize[0]
            chunk_x0 = icol*chunksize[1]
            chunk_x1 = (icol+1)*chunksize[1]

            # Compute region within chunk to place data in
            cy0 = max(yoff, yst) - chunk_y0
            cy1 = min(yoff+size[0], yen) - yst
            cx0 = max(xoff, xst) - chunk_x0
            cx1 = min(xoff+size[1], xen) - xst

            # Compute region to slice from data
            dy0 = max(0, yst-yoff)
            dy1 = min(yoff+size[0], chunk_y1)
            dx0 = max(0, xst-xoff)
            dx1 = min(xoff+size[1], chunk_x1)

            chunkdata[cy0:cy1, cx0:cx1] = array[dy0:dy1, dx0:dx1]

            # Return to data store
            self._store(chunkdata, i)
        return

    def getdata(self, xoff, yoff, size):
        result = np.empty(size, self.dtype)
        for i, yst, yen, xst, xen in self.getchunks(yoff, xoff, *size):

            # Compute the extents from the chunk to retain
            cy0 = max(yoff, yst) - yst
            cy1 = min(yoff+size[0], yen) - yst
            cx0 = max(xoff, xst) - xst
            cx1 = min(xoff+size[1], xen) - xst

            # Compute the bounds in the output
            oy0 = max(0, yst-yoff)
            oy1 = min(size[0], yen-yoff)
            ox0 = max(0, xst-xoff)
            ox1 = min(size[1], xen-xoff)

            result[oy0:oy1, ox0:ox1] = self._retrieve(i)[cy0:cy1, cx0:cx1]

        return result

