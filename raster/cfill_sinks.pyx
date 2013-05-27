""" Function to iteratively fill depressions in a raster surface. """

import numpy as np
cimport numpy as np
cimport cython

#@cython.profile(True)
@cython.boundscheck(False)
@cython.wraparound(False)
cdef int find_min(double[:] Z, int[:] FLAG):
    """ Return index of smallest item in Z where FLAG == 1 """
    cdef int itr, idx
    cdef double mn
    mn = 1e16
    for itr in range(len(Z)):
        if FLAG[itr] == 1:
            if Z[itr] < mn:
                idx = itr
                mn = Z[itr]
    return idx

#@cython.profile(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.always_allow_keywords(False)
cdef double[:] flatten(double[:,:] A):
    cdef double[:] Aflat
    cdef int i, j, idx
    Aflat = np.empty(A.size)
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            idx = i * A.shape[1] + j
            Aflat[idx] = A[i,j]
    return Aflat

#@cython.profile(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def fill_sinks_cy2(double[:,:] Z):
    """ Fill sinks in a DEM following the algorithm of Wang and Liu
    (2006).

        *Z*    :    2d array of elevation or potential data
                    (must not contain NaN!)

    Wang, L. and Liu, H. An efficient method for identifying and filling
    surface depressions in digital elevation models for hydrologic
    analysis and modelling. International Journal of Geographical
    Information Science, 20:2 (2006).
    """

    cdef int[:] OPENROW, OPENCOL, OPENIND, BROW, BCOL, CLOSED
    cdef double[:] ZFLAT, OPENZ, SPILL
    cdef int n, nx, ny, nopen, i, j      # counters and row/column indices
    cdef int bidx, oidx, nidx     # boundary, open, and neighbour indices
    cdef double Zn
    cdef tuple neighbours
    
    # Initialize SPILL, OPEN and CLOSED
    n = Z.shape[0] * Z.shape[1]
    ZFLAT = flatten(Z)
    SPILL = np.zeros(n, dtype=np.double)
    CLOSED = np.zeros(n, dtype=np.int32)
    OPENZ = np.zeros(n, dtype=np.double)
    OPENIND = np.zeros(n, dtype=np.int32)
    BROW = np.zeros(n, dtype=np.int32)
    BCOL = np.zeros(n, dtype=np.int32)

    # Get the boundary cells
    ny, nx = (Z.shape[0], Z.shape[1])
    bidx = 0
    for i in range(ny):
        BROW[bidx] = i
        BCOL[bidx] = 0
        bidx += 1
        BROW[bidx] = i
        BCOL[bidx] = nx-1
        bidx += 1

    for j in range(nx):
        BROW[bidx] = 0
        BCOL[bidx] = j
        bidx += 1
        BROW[bidx] = ny-1
        BCOL[bidx] = j
        bidx += 1

    # Add boundary to OPEN
    nopen = 0
    for i in range(bidx):
        oidx = BROW[i]*nx + BCOL[i]
        OPENZ[oidx] = ZFLAT[oidx]
        OPENIND[oidx] = 1
        nopen += 1
        SPILL[oidx] = ZFLAT[oidx]
    
    # Iteratively fill sinks
    while nopen > 0:

        oidx = find_min(OPENZ, OPENIND)
        OPENIND[oidx] = 0
        nopen -= 1
        CLOSED[oidx] = 1

        neighbours = (oidx-nx-1, oidx-nx, oidx-nx+1,
                      oidx-1,             oidx+1,
                      oidx+nx-1, oidx+nx, oidx+nx+1)

        for nidx in neighbours:

            if 0 <= nidx < n:
                if CLOSED[nidx] == 0 and OPENIND[nidx] == 0:
    
                    Zn = ZFLAT[nidx]
                    SPILL[nidx] = max(Zn, SPILL[oidx])
    
                    OPENZ[nidx] = Zn
                    OPENIND[nidx] = 1
                    nopen += 1

    return np.reshape(SPILL, (ny, nx))
