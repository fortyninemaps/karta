from math import ceil, sqrt
import numpy as np
cimport numpy as np
cimport cython

DTYPE_float64 = np.float64
ctypedef np.float64_t DTYPE_float64_t

DTYPE_int32 = np.int32
ctypedef np.int32_t DTYPE_int32_t

@cython.wraparound(False)
def fillarray_double(double[:,:] array not None,
                     int[:] I not None,
                     int[:] J not None,
                     double[:] Z not None,
                     double nodata_value):
    """ fillarray_double places values from *Z* into 2d *array* at the indices
    *I*, *J*. *array* should be a zero array.

    Returns 0 on success and 1 on failure.
    """
    cdef int ny, nx
    cdef int i, j, idx
    cdef double z
    cdef np.ndarray[DTYPE_int32_t, ndim=2] counts

    if (len(I) != len(J)) or (len(I) != len(Z)):
        return 1

    ny = array.shape[0]
    nx = array.shape[1]
    counts = np.zeros([ny, nx], dtype=DTYPE_int32)

    for idx in range(len(I)):
        i = I[idx]
        j = J[idx]
        array[i,j] += Z[idx]
        counts[i,j] += 1

    for i in range(ny):
        for j in range(nx):
            if counts[i,j] != 0:
                array[i,j] = array[i,j] / counts[i,j]
            else:
                array[i,j] = nodata_value
    return 0

cdef float interpolate1(float x, float y, float a, float b, float c, float d):
    """ Return a value *v(x,y)* in the regular structured stencil

            a --- b
            |  v  |
            c --- d

    using linear interpolation. The coordinates (x, y) must be normalized by
    the horizontal and vertical grid spacings, respectively.
    """
    cdef float left, right, res

    left = (c-a)*y + a
    right = (d-b)*y + b
    res = (right - left) * x + left
    return res

