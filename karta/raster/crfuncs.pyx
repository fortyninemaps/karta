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
    
def streamline2d(np.ndarray U, np.ndarray V, float x0, float y0,
                 float ds=0.5, int max_nodes=5000, tuple res=(1.0, 1.0),
                 float tol=0.0, int momentum=False):
    cdef int m, n, i
    cdef list X, Y
    cdef float dxds, dyds
    cdef float k1x, k1y, k2x, k2y, k3x, k3y, k4x, k4y
    cdef tuple result
    
    m = U.shape[0]
    n = U.shape[1]

    # Initialize streamline vectors
    X = [x0]
    Y = [y0]

    # Adjust the input arrays to make them have resolution 1x1
    U = U / res[0]
    V = V / res[1]

    # Put x0, y0 into grid units
    x0 = x0 / res[0]
    y0 = y0 / res[1]

    # Check that x0, y0 are within the bounds of U, V
    if (x0 >= n-1) or (y0 >= m-1) or (x0 < 0) or (y0 < 0):
        raise IndexError("starting position is outside vector field")

    i = 0

    while i < max_nodes:

        try:

            # Use linear interpolation to find dx/ds and dy/ds at (x0, y0)
            dxds = interpolate1(x0 % 1.0, y0 % 1.0,
                        U[y0, x0], U[y0, ceil(x0)],
                        U[ceil(y0), x0], U[ceil(y0), ceil(x0)])
            dyds = interpolate1(x0 % 1.0, y0 % 1.0,
                        V[y0, x0], V[y0, ceil(x0)],
                        V[ceil(y0), x0], V[ceil(y0), ceil(x0)])

            k1x = ds * dxds
            k1y = ds * dyds

            dxds = interpolate1(x0 % 1.0 + 0.5*ds, y0 % 1.0 + 0.5*k1x,
                        U[y0, x0], U[y0, ceil(x0)],
                        U[ceil(y0), x0], U[ceil(y0), ceil(x0)])
            dyds = interpolate1(x0 % 1.0 + 0.5*ds, y0 % 1.0 + 0.5*k1y,
                        V[y0, x0], V[y0, ceil(x0)],
                        V[ceil(y0), x0], V[ceil(y0), ceil(x0)])

            k2x = ds * dxds
            k2y = ds * dyds

            dxds = interpolate1(x0 % 1.0 + 0.5*ds, y0 % 1.0 + 0.5*k2x,
                        U[y0, x0], U[y0, ceil(x0)],
                        U[ceil(y0), x0], U[ceil(y0), ceil(x0)])
            dyds = interpolate1(x0 % 1.0 + 0.5*ds, y0 % 1.0 + 0.5*k2y,
                        V[y0, x0], V[y0, ceil(x0)],
                        V[ceil(y0), x0], V[ceil(y0), ceil(x0)])

            k3x = ds * dxds
            k3y = ds * dyds

            dxds = interpolate1(x0 % 1.0 + ds, y0 % 1.0 + k3x,
                        U[y0, x0], U[y0, ceil(x0)],
                        U[ceil(y0), x0], U[ceil(y0), ceil(x0)])
            dyds = interpolate1(x0 % 1.0 + ds, y0 % 1.0 + k3y,
                        V[y0, x0], V[y0, ceil(x0)],
                        V[ceil(y0), x0], V[ceil(y0), ceil(x0)])

            k4x = ds * dxds
            k4y = ds * dyds

            x0 = x0 + k1x/6.0 + k2x/3.0 + k3x/3.0 + k4x/6.0
            y0 = y0 + k1y/6.0 + k2y/3.0 + k3y/3.0 + k4y/6.0

        except IndexError:      # x0 or y0 was an invalid (e.g. nan)
            break

        # Check that x0, y0 are within the bounds of U, V
        if (x0 >= n-1) or (y0 >= m-1) or (x0 < 0) or (y0 < 0):
            break

        X.append(x0*res[0])
        Y.append(y0*res[1])

        # Check that rate of change is greater than the tolerance limit
        if tol > 0.0:
            if sqrt( k4x**2 + k4y**2 ) < tol:
                break

        i += 1

    result = (X, Y)
    return result
    
    
