"""
Flow-routing algorithms


The following license applies to:
    pixel_flow, facet_flow

    Copyright (c) 2009, The MathWorks, Inc.
    All rights reserved.

    Redistribution and use in source and binary forms, with or
    without modification, are permitted provided that the following
    conditions are met:

        * Redistributions of source code must retain the above copyright
          notice, this list of conditions and the following disclaimer.

        * Redistributions in binary form must reproduce the above copyright
          notice, this list of conditions and the following
          disclaimer in the documentation and/or other materials
          provided with the distribution

        * Neither the name of the The MathWorks, Inc. nor the names
          of its contributors may be used to endorse or promote
          products derived from this software without specific prior
          written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
    CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
    INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
    CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
    USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
    LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
    IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
    THE POSSIBILITY OF SUCH DAMAGE.

"""

import math
import numpy as np
from . import raster
from scipy import sparse

try:
    from . import crfuncs as rfuncs
except ImportError:
    from . import rfuncs as rfuncs

def facet_flow(e0, e1, e2, d1=1.0, d2=1.0):
    """ Return flow direction and slope for an east-northeast grid
    facet using the method of Tarboton (1997).

                        e2
                        |
                        |
                        | d2
                        |
                        |
            e0 -------- e1
                  d1

    Inspired by blog post by Steven L. Eddins (2007)
    """

    s1 = (e0 - e1) / d1             # Equation (1)
    s2 = (e1 - e2) / d2             # Equation (2)

    r = math.atan2(s2, s1)          # Equation (3)
    s = math.sqrt(s1*s1 + s2*s2)    # Equation (3)

    # Constrain flow direction to lie within or along the edges of the facet
    r, s = r > 0.0 and (r, s) or (0.0, s1)

    diagonal_angle = math.atan2(d2, d1)
    diagonal_distance = math.sqrt(d1*d1 + d2*d2)
    r, s = r < diagonal_angle and \
            (r, s) or (diagonal_angle, (e0 - e2) / diagonal_distance)

    return r, s


def pixel_flow(E, i, j, d1=1.0, d2=1.0):
    """ Downslope flow direction for DEM pixels using method of Tarboton (1997).

    Parameters
    ----------

    E : elevation array

    i : pixel row

    j : pixel column

    d1 : row spacing (default 1.0)

    d2 : column spacing (default 1.0)

    Inspired by blog post by Steven L. Eddins (2007)
    """

    def border_nans(A):
        """ Return indices of nans that touch the border. """
        #nans = np.isnan(A)
        raise Exception("border_nans not implemented")

    m, n = E.shape

    # Compute linear indices at desired locations
    e0_idx = i*(n-0) + j

    # Table 1, page 311
    # Row and column offsets corresponding to e1 and e2 for each entry
    e1_row_offsets = np.array([0, -1, -1,  0,  0,  1,  1,  0])
    e1_col_offsets = np.array([1,  0,  0, -1, -1,  0,  0,  1])

    e2_row_offsets = np.array([-1, -1, -1, -1,  1,  1,  1,  1])
    e2_col_offsets = np.array([ 1,  1, -1, -1, -1, -1,  1,  1])

    # Linear e1 and e2 offsets
    e1_linear_offsets = e1_row_offsets * (n-0) + e1_col_offsets
    e2_linear_offsets = e2_row_offsets * (n-0) + e2_col_offsets

    # Initialize R and S based on the first facet
    Eflat = E.flatten()
    E0 = Eflat[e0_idx]

    E1 = Eflat[e0_idx + e1_linear_offsets[0]]
    E2 = Eflat[e0_idx + e2_linear_offsets[0]]

    R, S = facet_flow(E0, E1, E2, d1, d2)

    # Multipliers ac and af used to convert from east-northeast facet
    # to any other fact (using Eq. 6)
    ac = [0,  1,  1,  2,  2,  3,  3,  4]
    af = [1, -1,  1, -1,  1, -1,  1, -1]

    R = (af[0] * R) + (ac[0] * math.pi / 2.0) if S > 0 else np.nan

    # Compute Rk and Sk corresponding to the k-th facet
    # Where Sk>0 and Sk>S, S=Sk and recompute R(Rk)
    for k in range(1,8):
        E1 = Eflat[e0_idx + e1_linear_offsets[k]]
        E2 = Eflat[e0_idx + e2_linear_offsets[k]]
        Rk, Sk = facet_flow(E0, E1, E2, d1, d2)

        if (Sk > S) and (Sk > 0.0):
            R = (af[k] * Rk) + (ac[k] * math.pi / 2.0)
        S = max(S, Sk)

    return R, S


def dem_flow(D):
    """ Calculate a flow field (aspect and slope) for an array *D*. Return
    aspect (*R*) and slope (*S*).

    Uses the D-infinity method of Tarboton (1997)
    Based on code implemented by Steven L. Eddins (2007)
    """
    Dp = raster.pad(D)
    c = [(i,j) for i in range(1, Dp.shape[0]-1) for j in range(1, Dp.shape[1]-1)]

    ff = map(lambda a: pixel_flow(Dp, a[0], a[1], d1=1.0, d2=1.0), c)
    R = np.array([i[0] for i in ff]).reshape(D.shape)
    S = np.array([i[1] for i in ff]).reshape(D.shape)
    return R, S

def dem_flow2(D):
    """ Calculate aspect and slope using karta.raster functions """
    pi = np.pi
    asp = raster.aspect(D)
    asp[0,:] = 0.5*pi
    asp[-1,:] = 1.5*pi
    asp[:,0] = pi
    asp[:,-1] = 0.0
    asp[0,0] = 0.75*pi
    asp[0,-1] = 0.25*pi
    asp[-1,0] = 1.25*pi
    asp[-1,-1] = 1.75*pi
    slope = raster.slope(D)
    return asp, slope

def prop_dinfty(position, alpha):
    """ Use the D-oo algorithm to return proportion of flow to the center of a
    3x3 grid.

    Parameters:
    -----------

    position : cell location assigned based on the map

                    0   1   2
                    3   *   4
                    5   6   7

    alpha : angle measured counter-clockwise from a vector pointing toward
    position 4
    """
    pi = math.pi

    # beta is the direction to the centre node
    BETA = [1.75, 1.5, 1.25, 0.0, 1.0, 0.25, 0.5, 0.75]
    beta = BETA[position] * pi

    theta = abs(rfuncs.diffrad(alpha, beta))
    if theta <= 0.25*pi and np.isnan(alpha) == False:
        P = (0.25*pi-theta) / (0.25*pi)
    else:
        P = 0.0
    return P


def prop_d8(position, alpha):
    """ Use the D8 algorthm to return proportion of flow to the center of a 3x3
    grid.

    Parameters:
    -----------

    position : cell location assigned based on the map

                    0   1   2
                    3   *   4
                    5   6   7

    alpha : angle measured counter-clockwise from a vector pointing toward
    position 4
    """
    pi = math.pi

    # beta is the direction to the centre node
    BETA = [1.75, 1.5, 1.25, 0.0, 1.0, 0.25, 0.5, 0.75]
    beta = BETA[position] * pi

    theta = abs(rfuncs.diffrad(alpha, beta))
    if theta <= 0.125*pi and np.isnan(alpha) == False:
        P = 1.0
    else:
        P = 0.0
    return P


def upslope_area(F, A, proportion=prop_dinfty):
    """ Calculate upslope area with a sparse matrix formulation.

    Parameters:
    -----------

    F : flow direction array (e.g. from dem_flow())

    A : array providing the area of each grid cell

    proportion : function that defines how flow should be partitioned (e.g. D8,
        D-infinity algorithms).
    """
    m, n = F.shape
    if not hasattr(A, 'shape'):
        A = np.ones_like(F) * A

    # Assemble a flow contribution matrix
    rows = np.zeros(F.size*9)
    cols = np.zeros(F.size*9)
    data = np.zeros(F.size*9)

    cnt = 0
    for i in range(m):
        for j in range(n):
            i_ = i*n + j
            i0 = i_ - n - 1
            i1 = i_ - n
            i2 = i_ - n + 1
            i3 = i_ - 1
            i4 = i_ + 1
            i5 = i_ + n - 1
            i6 = i_ + n
            i7 = i_ + n + 1

            rows[cnt] = i_
            cols[cnt] = i_
            data[cnt] = A[i,j]
            cnt += 1

            if i > 0:
                if j > 0:
                    rows[cnt] = i_
                    cols[cnt] = i0
                    data[cnt] = -A[i-1,j-1] * proportion(0, F[i-1,j-1])
                    cnt += 1

                rows[cnt] = i_
                cols[cnt] = i1
                data[cnt] = -A[i-1,j]   * proportion(1, F[i-1,j])
                cnt += 1

                if j < n-1:
                    rows[cnt] = i_
                    cols[cnt] = i2
                    data[cnt] = -A[i-1,j+1] * proportion(2, F[i-1,j+1])
                    cnt += 1

            if j > 0:
                rows[cnt] = i_
                cols[cnt] = i3
                data[cnt] = -A[i,j-1]   * proportion(3, F[i,j-1])
                cnt += 1

            if j < n-1:
                rows[cnt] = i_
                cols[cnt] = i4
                data[cnt] = -A[i,j+1]   * proportion(4, F[i,j+1])
                cnt += 1

            if i < m-1:
                if j > 0:
                    rows[cnt] = i_
                    cols[cnt] = i5
                    data[cnt] = -A[i+1,j-1] * proportion(5, F[i+1,j-1])
                    cnt += 1

                rows[cnt] = i_
                cols[cnt] = i6
                data[cnt] = -A[i+1,j]   * proportion(6, F[i+1,j])
                cnt += 1

                if j < n-1:
                    rows[cnt] = i_
                    cols[cnt] = i7
                    data[cnt] = -A[i+1,j+1] * proportion(7, F[i+1,j+1])
                    cnt += 1

    C = sparse.coo_matrix((data, (rows, cols)), shape=(F.size, F.size))

    # Assemble a flow contribution matrix
    #C = sparse.lil_matrix((F.size, F.size))
    #for i in xrange(m):
    #    for j in xrange(n):

    #        i_ = i*n + j
    #        i0 = i_ - n - 1
    #        i1 = i_ - n
    #        i2 = i_ - n + 1
    #        i3 = i_ - 1
    #        i4 = i_ + 1
    #        i5 = i_ + n - 1
    #        i6 = i_ + n
    #        i7 = i_ + n + 1

    #        C[i_,i_] = A[i,j]
    #        if i > 0:
    #            if j > 0:
    #                C[i_,i0] = -A[i-1,j-1] * proportion(0, F[i-1,j-1])
    #            C[i_,i1]     = -A[i-1,j]   * proportion(1, F[i-1,j])
    #            if j < n-1:
    #                C[i_,i2] = -A[i-1,j+1] * proportion(2, F[i-1,j+1])
    #        if j > 0:
    #            C[i_,i3]     = -A[i,j-1]   * proportion(3, F[i,j-1])
    #        if j < n-1:
    #            C[i_,i4]     = -A[i,j+1]   * proportion(4, F[i,j+1])
    #        if i < m-1:
    #            if j > 0:
    #                C[i_,i5] = -A[i+1,j-1] * proportion(5, F[i+1,j-1])
    #            C[i_,i6]     = -A[i+1,j]   * proportion(6, F[i+1,j])
    #            if j < n-1:
    #                C[i_,i7] = -A[i+1,j+1] * proportion(7, F[i+1,j+1])

    # Next, solve the linear problem: C * U_a = \summation{A}
    U = sparse.linalg.spsolve(C.tocsc(), A.flatten() * np.ones(m*n))
    return U.reshape(F.shape)

#def upslope_area_nn(F, A, proportion=prop_dinfty):
#    """ Calculate upslope area with a sparse matrix formulation.
#
#    Parameters:
#    -----------
#
#    F : flow direction array (e.g. from dem_flow())
#
#    A : array providing the area of each grid cell
#
#    proportion : function that defines how flow should be partitioned (e.g. D8,
#        D-infinity algorithms).
#    """
#    m, n = F.shape
#    if not hasattr(A, 'shape'):
#        A = np.ones_like(F) * A
#
#    nnmap = np.nonzero(np.isnan(F) == False)
#    nnsize = len(nnmap[0])
#
#    # Assemble a flow contribution matrix
#    C = sparse.lil_matrix((nnsize+2, nnsize+2))
#    for cell in range(nnsize):
#        i = nnmap[0][cell] + 1
#        j = nnmap[1][cell] + 1
#
#        i_ = i*n + j
#        i0 = i_ - n - 1
#        i1 = i_ - n
#        i2 = i_ - n + 1
#        i3 = i_ - 1
#        i4 = i_ + 1
#        i5 = i_ + n - 1
#        i6 = i_ + n
#        i7 = i_ + n + 1
#
#        C[i_,i_] = A[i,j]
#        if i > 0:
#            if j > 0:
#                C[i_,i0] = -A[i-1,j-1] * proportion(0, F[i-1,j-1])
#            C[i_,i1]     = -A[i-1,j]   * proportion(1, F[i-1,j])
#            if j < n-1:
#                C[i_,i2] = -A[i-1,j+1] * proportion(2, F[i-1,j+1])
#        if j > 0:
#            C[i_,i3]     = -A[i,j-1]   * proportion(3, F[i,j-1])
#        if j < n-1:
#            C[i_,i4]     = -A[i,j+1]   * proportion(4, F[i,j+1])
#        if i < m-1:
#            if j > 0:
#                C[i_,i5] = -A[i+1,j-1] * proportion(5, F[i+1,j-1])
#            C[i_,i6]     = -A[i+1,j]   * proportion(6, F[i+1,j])
#            if j < n-1:
#                C[i_,i7] = -A[i+1,j+1] * proportion(7, F[i+1,j+1])
#
#    # Next, solve the linear problem: C * U_a = \summation{A}
#    Unn = sparse.linalg.spsolve(C.tocsc(), A.flatten())
#
#    # Reconstruct the full array
#    U = np.zeros_like(F)
#    Unndiag = Unn.diagonal()
#    for cell in range(nnsize):
#        i = nnmap[0][cell] + 1
#        j = nnmap[1][cell] + 1
#        U[i,j] = Unnfull[cell]
#
#    return U

