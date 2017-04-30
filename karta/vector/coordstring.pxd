cimport numpy as np
from cpython cimport bool

cdef class CoordString:
    cdef int length
    cdef double *coords
    cdef readonly int rank
    cdef readonly bool ring
    cpdef np.ndarray slice(self, int start, int stop=?, int step=?)
    cdef double getX(self, int index)
    cdef double getY(self, int index)
    cdef double getZ(self, int index)
