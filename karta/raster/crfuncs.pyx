from math import ceil, sqrt
import numpy as np
cimport numpy as np
from cpython cimport bool

def diffrad(double a, double b):
    """ Return the difference in radians between two angles """
    cdef double pi
    pi = 3.141592653589793
    return (a + pi - b) % (2*pi) - pi

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
    
    
