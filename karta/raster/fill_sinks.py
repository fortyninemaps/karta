""" Function to iteratively fill depressions in a raster surface. """

import heapq
import numpy as np

def neighbours_of(a):
    """ For a (z,i,j) point `a`, return the neighbouring indices. """
    _, i, j = a
    return ((i-1, j-1), (i, j-1), (i+1, j-1),
            (i-1, j), (i+1, j),
            (i-1, j+1), (i, j+1), (i+1, j+1))

def fill_sinks(Z):
    """ Fill sinks in a DEM following the algorithm of Wang and Liu
    (2006).

        *Z*    :    2d array of elevation or potential data
                    (must not contain NaN!)

    Wang, L. and Liu, H. An efficient method for identifying and filling
    surface depressions in digital elevation models for hydrologic
    analysis and modelling. International Journal of Geographical
    Information Science, 20:2 (2006).
    """

    # Initialize SPILL and CLOSED
    SPILL = Z.copy()
    CLOSED = np.zeros_like(Z)
    OPEN = []
    nopen = 0

    # Get the boundary cells
    ny, nx = (Z.shape[0], Z.shape[1])
    B = [(i, j) for i in range(ny) for j in (0, nx-1)]
    B.extend([(i, j) for i in (0, ny-1) for j in range(1, nx-1)])

    # Get z along the boundary
    for b in B:
        SPILL[b] = Z[b]
        heapq.heappush(OPEN, (Z[b], b[0], b[1]))
        nopen += 1

    while nopen > 0:

        c = heapq.heappop(OPEN)
        nopen -= 1
        CLOSED[c[1:]] = 1

        for n in neighbours_of(c):
            if (0 <= n[0] < ny) and (0 <= n[1] < nx):
                Zn = Z[n]
                if not CLOSED[n] and (Zn, n[0], n[1]) not in list(OPEN):
                    SPILL[n] = max(Zn, SPILL[c[1], c[2]])
                    heapq.heappush(OPEN, (Zn, n[0], n[1]))
                    nopen += 1

    return SPILL
