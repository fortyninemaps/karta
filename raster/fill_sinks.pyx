import numpy as np
from itertools import count, izip
cimport numpy as np

def neighbours_of(tuple a not None):
    cdef int i, j, k
    i, j, z = a
    return ((i-1, j-1), (i, j-1), (i+1, j-1),
            (i-1, j), (i+1, j),
            (i-1, j+1), (i, j+1), (i+1, j+1))

def fill_sinks(np.ndarray Z not None):
    """ Fill sinks in a DEM following the algorithm of Wang and Liu
    (2006).

        *Z*    :    2d array of elevation or potential data
                    (must not contain NaN!)

    Wang, L. and Liu, H. An efficient method for identifying and filling
    surface depressions in digital elevation models for hydrologic
    analysis and modelling. International Journal of Geographical
    Information Science, 20:2 (2006).
    """

    cdef np.ndarray SPILL, CLOSED
    cdef list OPEN, B
    cdef int nx, ny, i, j, a
    cdef tuple b, c
    cdef tuple N, n
    cdef float Zn, lowest
    
    
    # Initialize SPILL and CLOSED
    SPILL = Z.copy()
    CLOSED = np.zeros_like(Z)
    OPEN = []

    # Get the boundary cells
    ny, nx = (Z.shape[0], Z.shape[1])
    B = [(i, j) for i in range(ny) for j in (0, nx-1)]
    B.extend([(i, j) for i in (0, ny-1) for j in range(nx)])

    # Get z along the boundary
    for b in B:
        SPILL[b] = Z[b]
        OPEN.append((b[0], b[1], Z[b]))

    while len(OPEN) > 0:

        #OPEN.sort(key=lambda a: a[2], reverse=True)
        #c = OPEN.pop()
        # This fancy version uses iterators and generators
        c = OPEN.pop(min(izip((c[2] for c in OPEN), count()))[1])
        CLOSED[c[:2]] = 1

        for n in neighbours_of(c):
            #if (n[0]>=ny or n[0]<0 or n[1]>=nx or n[1]<0) is False:
            if (n[0]<ny) and (n[0]>=0) and (n[1]<nx) and (n[1]>=0):
                Zn = Z[n]
                if CLOSED[n] == 0 and (n[0], n[1], Zn) not in OPEN:
                    SPILL[n] = max(Zn, SPILL[c[0], c[1]])
                    OPEN.append((n[0], n[1], Zn))

    return SPILL
