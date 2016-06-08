import numpy as np
cimport numpy as np
from cpython cimport bool
from cpython.array cimport array, clone

cdef double mind(double a, double b):
    return a if a<= b else b

cdef double maxd(double a, double b):
    return a if a>= b else b

cdef int floord(double a):
    if (a % 1) == 0:
        return <int> a
    else:
        return <int> (a // 1)

cdef int ceild(double a):
    if (a % 1) == 0:
        return <int> a
    else:
        return floord(a) + 1

cdef array template_dbl = array("d", [])

cdef class CoordString:

    cdef readonly double[:] coords
    cdef readonly int rank

    def __cinit__(self, object coords):
        cdef int length = len(coords)
        if length == 0:
            self.rank = 0
        else:
            self.rank = len(coords[0])

        if self.rank not in (0, 2, 3):
            raise ValueError("non-empty Geometry rank must be 2 or 3")

        cdef array coords_array = clone(template_dbl, length*self.rank, False)
        self.coords = coords_array

        cdef int i, j
        cdef object xy
        i = 0
        for xy in coords:
            for j in range(self.rank):
                self.coords[i+j] = xy[j]
            i = i + self.rank
        return

    def __len__(self):
        if self.rank == 0:
            return 0
        else:
            return int(len(self.coords)/self.rank)

    def __iter__(self):
        cdef int i = 0, n = len(self)
        while i < n:
            yield self[i]
            i += 1

    def __getitem__(self, int index):
        cdef double x, y, z
        cdef int pos = index*self.rank
        x = self.coords[pos]
        y = self.coords[pos+1]
        if self.rank == 2:
            return x, y
        else:
            z = self.coords[pos+2]
            return x, y, z

    def __setitem__(self, int key, double[:] value):
        cdef int idx = key*self.rank
        self.coords[key*self.rank:(key+1)*self.rank] = value

    def __hash__(self):
        return hash(self.coords)

    def __richcmp__(self, other, int op):
        if op == 2:
            return self._ceq(other)

    cdef bool _ceq(self, CoordString other):
        """ Tests equality """
        cdef int len1 = len(self.coords), len2 = len(other.coords)
        cdef int i = 0
        if len1 == len2:
            for i in range(len1):
                if self.coords[i] != other.coords[i]:
                    return False
            return True
        return False

    cpdef np.ndarray slice(self, int start, int stop=0, int step=1):
        """ Slice coordinate string, returning an <n x rank> numpy array. """
        while start < 0:
            start += len(self)
        while stop <= 0:
            stop += len(self)

        cdef int outlength
        if step != 0:
            outlength = ceild(<double> (abs(stop) - abs(start)) / abs(step))
        else:
            raise ValueError("step cannot equal zero")

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

        if len(self) == 0:
            return (np.nan, np.nan, np.nan, np.nan)

        xmin = self.coords[0]
        xmax = self.coords[0]
        ymin = self.coords[1]
        ymax = self.coords[1]
        while i != len(self)-1:
            offset += self.rank
            xmin = mind(xmin, self.coords[offset])
            xmax = maxd(xmax, self.coords[offset])
            ymin = mind(ymin, self.coords[offset+1])
            ymax = maxd(ymax, self.coords[offset+1])
            i += 1
        return (xmin, ymin, xmax, ymax)

    def vectors(self):
        cdef int i = 0, pos = 0, n = len(self)
        cdef np.ndarray x, y, z
        x = np.empty(n, dtype=np.float64)
        y = np.empty(n, dtype=np.float64)
        while i < n:
            x[i] = self.coords[pos]
            y[i] = self.coords[pos+1]
            i += 1
            pos += self.rank
        if self.rank == 2:
            return x, y
        else:
            z = np.empty(n, dtype=np.float64)
            pos = 0
            i = 0
            while i < n:
                z[i] = self.coords[pos+2]
                i += 1
                pos += self.rank
            return x, y, z

    def asarray(self):
        return np.array(self.coords).reshape([-1, self.rank])
