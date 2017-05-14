from libc.stdlib cimport malloc, free
from libc.math cimport isnan, NAN
import numpy as np
cimport numpy as np
from cpython cimport bool

cdef double mind(double a, double b) nogil:
    return a if a<= b else b

cdef double maxd(double a, double b) nogil:
    return a if a>= b else b

cdef int floord(double a) nogil:
    if (a % 1) == 0:
        return <int> a
    else:
        return <int> (a // 1)

cdef int ceild(double a) nogil:
    if (a % 1) == 0:
        return <int> a
    else:
        return floord(a) + 1

cdef class CoordString:
    """ A CoordString is a datastructure for coordinate data. Initialize a
    CoordString with 2 or 3 iterable objects, corresponding to x, y, and
    optionally z coordinates.

    Parameters
    ----------
    coords : iterable
        list of coordinates, with items of length 2 or 3
    ring : bool
        indicates whether coordinate are implicitly closed
    """
    def __cinit__(self, object coords, bool ring=False):
        cdef int length = -1
        try:
            length = len(coords)
            if length == 0:
                self.rank = -1
            else:
                self.rank = len(coords[0])
        except TypeError:
            raise ValueError("coords argument must be an iterable of iterables")

        if self.rank not in (-1, 2, 3):
            raise ValueError("CoordString rank must be 2 or 3")

        self.ring = ring
        self.length = length
        self.coords = <double *> malloc(length * self.rank * sizeof(double))

        cdef int i = 0, j
        cdef object xy
        for xy in coords:
            for j in range(self.rank):
                if isnan(xy[j]):
                    raise ValueError("coordinate object contains NaN values")
                self.coords[i+j] = xy[j]
            i = i + self.rank
        return

    def __dealloc__(self):
        free(self.coords)

    def __len__(self):
        if self.rank == -1:
            return 0
        else:
            return self.length

    def __iter__(self):
        cdef int i = 0
        while i < self.length:
            yield self[i]
            i += 1

    def __getitem__(self, int index):
        if index > self.length:
            raise IndexError("index out of range")
        if self.rank == -1:
            raise IndexError("CoordString empty")

        cdef double x, y, z
        cdef int pos = (index*self.rank) % (self.length*self.rank)
        x = self.coords[pos]
        y = self.coords[pos+1]
        if self.rank == 2:
            return x, y
        else:
            z = self.coords[pos+2]
            return x, y, z

    def __setitem__(self, int key, double[:] value):
        cdef int idx = key*self.rank
        cdef int i = 0
        if self.rank == -1:
            raise IndexError("CoordString empty")
        for i in range(self.rank):
            self.coords[idx+i] = value[i]

    def __hash__(self):
        cdef np.ndarray[np.double_t, ndim=1] buf = np.empty(self.length, dtype=np.double)
        cdef int i = 0
        for i in range(self.length):
            buf[i] = self.coords[i]
        return hash(buf.tostring())

    def __richcmp__(self, other, int op):
        if op == 2:
            return self._ceq(other)
        elif op == 3:
            return not self._ceq(other)
        else:
            raise NotImplementedError("CoordString comparison")

    def _ceq(self, CoordString other):
        """ Tests equality """
        cdef int i = 0
        if self.length != other.length:
            return False
        for i in range(self.length):
            if self.coords[i] != other.coords[i]:
                return False
        return True

    cdef double getX(self, int index):
        if index == self.length and self.ring:
            index = 0
        return self.coords[index*self.rank]

    cdef double getY(self, int index):
        if index == self.length and self.ring:
            index = 0
        return self.coords[index*self.rank+1]

    cdef double getZ(self, int index):
        if self.rank != 3:
            return NAN
        if index == self.length and self.ring:
            index = 0
        return self.coords[index*self.rank+2]

    cpdef np.ndarray slice(self, int start, int stop=0, int step=1):
        """ Slice coordinate string, returning an <n x rank> numpy array. """
        if self.rank == -1:
            return ValueError("CoordString empty")

        while start < 0:
            start += self.length
        while stop <= 0:
            stop += self.length

        cdef int outlength
        if step != 0:
            outlength = ceild(<double> (abs(stop) - abs(start)) / abs(step))
        else:
            raise ValueError("step cannot equal zero")

        cdef np.ndarray[np.double_t, ndim=1] result = np.empty(outlength*self.rank, dtype=np.double)
        cdef int i = 0, j = 0, pos = start
        cdef int rank = self.rank
        while pos < stop:
            for j in range(self.rank):
                result[i] = self.coords[pos*rank+j]
                i += 1
            pos += step
        return result.reshape([outlength, self.rank])

    @property
    def bbox(self):
        cdef double xmin, ymin, xmax, ymax
        cdef int i = 0
        cdef int offset = 0

        if self.length == 0:
            return (NAN, NAN, NAN, NAN)

        xmin = self.coords[0]
        xmax = self.coords[0]
        ymin = self.coords[1]
        ymax = self.coords[1]
        while i != self.length-1:
            offset += self.rank
            xmin = mind(xmin, self.coords[offset])
            xmax = maxd(xmax, self.coords[offset])
            ymin = mind(ymin, self.coords[offset+1])
            ymax = maxd(ymax, self.coords[offset+1])
            i += 1
        return (xmin, ymin, xmax, ymax)

    def vectors(self, bool drop_z=False):
        cdef int i = 0, pos = 0
        cdef np.ndarray[np.double_t, ndim=1] x, y
        x = np.empty(self.length, dtype=np.double)
        y = np.empty(self.length, dtype=np.double)
        while i < self.length:
            x[i] = self.coords[pos]
            y[i] = self.coords[pos+1]
            i += 1
            pos += self.rank
        if (self.rank == 2) or drop_z:
            return x, y

        cdef np.ndarray[np.double_t, ndim=1] z
        z = np.empty(self.length, dtype=np.double)
        pos = 0
        i = 0
        while i < self.length:
            z[i] = self.coords[pos+2]
            i += 1
            pos += self.rank
        return x, y, z

    def asarray(self):
        if self.rank == -1:
            return np.array([[]], dtype=np.double)
        cdef np.ndarray[np.double_t, ndim=1] arr = np.empty(self.rank*self.length, dtype=np.double)
        cdef int i = 0
        for i in range(self.length*self.rank):
            arr[i] = self.coords[i]
        return arr.reshape([-1, self.rank])
