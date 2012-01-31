"""
Submodule for generating streamlines from vector field rasters.
"""

from math import ceil

def streamline2d(U, V, x0, y0, ds=0.5, max_nodes=1000, res=(1.0, 1.0)):
    """ Integrate velocity field (*U*, *V*) using a 4th-order Runge-Kutte
    scheme, starting from *x0*, *y0*. A pair of lists with coordinates
    in X and Y is returned.

    *ds* is the step size
    *max_nodes* is the maximum number of steps to take. Iteration will
        be terminated automatically if the streamline reaches the edge
        of the vector field.
    *res* is the resolution of *U*, *V*, in the same units as *x0*, *y0*.
    """

    def interpolate1(x, y, a, b, c, d):
        left = (c-a)*y + a
        right = (d-b)*y + b
        return (right - left) * x + left

    m, n = U.shape

    # Initialize streamline vectors
    X = [x0]
    Y = [y0]

    # Put x0, y0 into grid units
    x0 = x0 / res[0]
    y0 = y0 / res[1]

    # Check that x0, y0 are within the bounds of U, V
    if (x0 >= n-1) or (y0 >= m-1) or (x0 < 0) or (y0 < 0):
        raise IndexError("starting position is outside vector field")

    i = 0

    while i < max_nodes:

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

        # Check that x0, y0 are within the bounds of U, V
        if (x0 >= n-1) or (y0 >= m-1) or (x0 < 0) or (y0 < 0):
            break

        X.append(x0*res[0])
        Y.append(y0*res[1])

        i += 1

    return X, Y


