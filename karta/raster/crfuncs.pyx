from math import ceil, sqrt
import numpy as np
cimport numpy as np
cimport cython

DTYPE_float64 = np.float64
ctypedef np.float64_t DTYPE_float64_t

DTYPE_int32 = np.int32
ctypedef np.int32_t DTYPE_int32_t

def get_positions_vec(tuple T, double[:] x not None, double[:] y not None):
    """ Return the row and column indices according to grid transform T for
    points in vectors x and y. """
    # (a, b, c, d, e, f) = T
    # X = a + jc + ie
    # Y = b + id + jf

    # i = (X - a - jc) / e
    # i = (Y - b - jf) / d
    # e(Y - b - jf) = d(X - a - jc)
    # ey - eb - ejf = dx - da - djc
    # djc - ejf = dx - da + eb - ey
    # j ( dc - ef ) = dx - da + eb - ey
    # j = (dx - da + eb - ey) / (dc - ef)
    cdef int i = 0
    cdef int n = len(x)
    cdef double x0, y0, dx, dy, sx, sy
    cdef double i_, j_
    cdef np.ndarray[np.float64_t] I = np.empty(n)
    cdef np.ndarray[np.float64_t] J = np.empty(n)
    x0 = T[0]
    y0 = T[1]
    dx = T[2]
    dy = T[3]
    sx = T[4]
    sy = T[5]
    while i != n:
        j_ = (dy*x[i] - dy*x0 + sx*y0 - sx*y[i]) / (dx*dy - sx*sy)
        i_ = (y[i] - y0 - j_*sy) / dy
        J[i] = j_ - 0.5
        I[i] = i_ - 0.5
        i += 1
    return I, J

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

