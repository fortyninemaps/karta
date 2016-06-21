cimport numpy as np

cdef class CoordString:
    cdef readonly double[:] coords
    cdef readonly int rank
    cpdef np.ndarray slice(self, int start, int stop=?, int step=?)
