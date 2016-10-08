import numpy as np
cimport numpy as np
from cpython cimport bool
from cpython.array cimport array, clone

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

cdef array template_dbl = array("d", [])

cdef class CoordString:
    """ A CoordString is a datastructure for coordinate data. Initialize a
    CoordString with 2 or 3 iterable objects, corresponding to x, y, and
    optionally z coordinates.

    Parameters
    ----------
    coords : iterable
        list of coordinates, with items of length 2 or 3
    """
    def __cinit__(self, object coords, bool ring=False):
        cdef int length = -1
        try:
            length = len(coords)
            self.length = length
            if length == 0:
                self.rank = -1
            else:
                self.rank = len(coords[0])
        except TypeError:
            raise ValueError("coords argument must be an iterable of iterables")

        if self.rank not in (-1, 2, 3):
            raise ValueError("CoordString rank must be 2 or 3")

        self.ring = ring

        cdef array coords_array = clone(template_dbl, length*self.rank, False)
        self.coords = coords_array

        cdef int i = 0, j
        cdef object xy
        for xy in coords:
            for j in range(self.rank):
                self.coords[i+j] = xy[j]
            i = i + self.rank
        return

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
        cdef double x, y, z
        cdef int pos = index*self.rank
        if self.rank == -1:
            raise IndexError("CoordString empty")
        x = self.coords[pos]
        y = self.coords[pos+1]
        if self.rank == 2:
            return x, y
        else:
            z = self.coords[pos+2]
            return x, y, z

    def __setitem__(self, int key, double[:] value):
        cdef int idx = key*self.rank
        if self.rank == -1:
            raise IndexError("CoordString empty")
        self.coords[key*self.rank:(key+1)*self.rank] = value

    def __hash__(self):
        return hash(self.coords)

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
        if self.ring and index == self.length:
            index = 0
        return self.coords[index*self.rank]

    cdef double getY(self, int index):
        if self.ring and index == self.length:
            index = 0
        return self.coords[index*self.rank+1]

    cdef double getZ(self, int index):
        if self.rank != 3:
            return np.nan
        if self.ring and index == self.length:
            index = 0
        return self.coords[index*self.rank+2]

    cpdef np.ndarray slice(self, int start, int stop=0, int step=1):
        """ Slice coordinate string, returning an <n x rank> numpy array. """
        while start < 0:
            start += self.length
        while stop <= 0:
            stop += self.length

        cdef int outlength
        if step != 0:
            outlength = ceild(<double> (abs(stop) - abs(start)) / abs(step))
        else:
            raise ValueError("step cannot equal zero")

        if self.rank == -1:
            return ValueError("CoordString empty")

        cdef np.ndarray result = np.empty([outlength, self.rank], dtype=np.float64)
        cdef int i = 0, pos = start
        cdef int rank = self.rank
        while pos < stop:
            result[i,:] = self.coords[pos*rank:(pos+1)*rank]
            i += 1
            pos += step
        return result

    @property
    def bbox(self):
        cdef double xmin, ymin, xmax, ymax
        cdef int i = 0
        cdef int offset = 0

        if self.length == 0:
            return (np.nan, np.nan, np.nan, np.nan)

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
        cdef np.ndarray x, y
        x = np.empty(self.length, dtype=np.float64)
        y = np.empty(self.length, dtype=np.float64)
        while i < self.length:
            x[i] = self.coords[pos]
            y[i] = self.coords[pos+1]
            i += 1
            pos += self.rank
        if (self.rank == 2) or drop_z:
            return x, y

        cdef np.ndarray z
        z = np.empty(self.length, dtype=np.float64)
        pos = 0
        i = 0
        while i < self.length:
            z[i] = self.coords[pos+2]
            i += 1
            pos += self.rank
        return x, y, z

    def asarray(self):
        if self.rank == -1:
            return np.array([[]], dtype=np.float64)
        return np.array(self.coords).reshape([-1, self.rank])
