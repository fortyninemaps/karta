"""
Flow-routing algorithms


The following license applies to:
    pixel_flow, facet_flow, upslope_area

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

#import sys
import raster
import math
import numpy as np
from scipy import sparse
import pdb, time


def facet_flow(e0, e1, e2, d1=30.0, d2=30.0):
    """ Return flow direction and slope for an east-northeast grid
    facet.

                        e2
                        |
                        |
                        | d2
                        |
                        |
            e0 -------- e1
                  d1

    Using method of Tarboton (1997)
    Based on code implemented by Steven L. Eddins (2007)
    """

    s1 = (e0 - e1) / d1             # Equation (1)
    s2 = (e1 - e2) / d2             # Equation (2)

    r = math.atan2(s2, s1)          # Equation (3)
    s = math.sqrt(s1*s1 + s2*s2)    # Equation (3)

    # Constrain flow direction to lie within or along the edges of the facet
    r, s = r > 0.0 and (r, s) or (0.0, s1)

    diagonal_angle = math.atan2(d2, d1)
    diagonal_distance = math.sqrt(d1*d1 + d2*d2)
    r, s = r < diagonal_angle and (r, s) or (diagonal_angle, (e0 - e2) / diagonal_distance)

    return r, s


def pixel_flow(E, i, j, d1=30.0, d2=30.0):
    """ Downslope flow direction for DEM pixels.

    Using method of Tarboton (1997)
    Based on code implemented by Steven L. Eddins (2007)
    """

    def border_nans(A):
        """ Return indices of nans that touch the border. """
        raise Exception("border_nans not implemented")
        nans = np.isnan(A)

    m, n = E.shape

    # Preprocess NaNs connected to the border.
    #highest_value = E[np.isnan(E)==False].max()
    #bump = 1.0
    #E[border_nans(E)] = highest_value + bump

    # Compute linear indices at desired locations
    e0_idx = j * m + i

    # Table 1, page 311
    # Row and column offsets corresponding to e1 and e2 for each entry
    e1_row_offsets = np.array([0, -1, -1,  0,  0,  1,  1,  0])
    e1_col_offsets = np.array([1,  0,  0, -1, -1,  0,  0,  1])

    e2_row_offsets = np.array([-1, -1, -1, -1,  1,  1,  1,  1])
    e2_col_offsets = np.array([ 1,  1, -1, -1, -1, -1,  1,  1])

    # Linear e1 and e2 offsets
    e1_linear_offsets = e1_col_offsets * m + e1_row_offsets
    e2_linear_offsets = e2_col_offsets * m + e2_row_offsets

    # Initialize R and S based on the first facet
    Eflat = E.flat
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
    """ Calculate a flow field (aspect and slope) for an array.

    Uses the D-infinity method of Tarboton (1997)
    Based on code implemented by Steven L. Eddins (2007)
    """
    Dp = raster.pad(D)
    c = [(i,j) for i in range(1, Dp.shape[0]-1) for j in range(1, Dp.shape[1]-1)]

    ff = map(lambda a: pixel_flow(Dp, a[0], a[1], d1=20.0, d2=20.0), c)
    R = np.array([i[0] for i in ff]).reshape(D.shape)
    S = np.array([i[1] for i in ff]).reshape(D.shape)

    return R, S


def build_contribution_matrix(R):
    """ Constructs an (n x n) sparse matrix P, where n is the number of
    elements in flow direction matrix R.
    """
    n = R.size
    nrow, ncol = R.shape
    pi = math.pi
    Rf = R.flat

    P = sparse.lil_matrix((n,n))

    # The primary diagonal
    P.setdiag(np.ones(n), k=0)

    for i in range(1,nrow-1):
        for j in range(1,ncol-1):
            m = i*nrow + j
            alpha = R[i,j]
            try:
                if alpha > 1.75*pi or alpha < 0.25*pi:
                    # Northward
                    P[m-ncol, m-2*ncol] = abs(0.25*pi - alpha) / 0.25*pi
                elif alpha > 0.25*pi and alpha < 0.75*pi:
                    # Eastward
                    P[m+1, m+2] = abs(0.25*pi - alpha + 0.5*pi) / 0.25*pi
                elif alpha > 0.75*pi and alpha < 1.25*pi:
                    # Southward
                    P[m+ncol, m+2*ncol] = abs(0.25*pi - alpha + pi) / 0.25*pi
                elif alpha > 1.25*pi and alpha < 1.75*pi:
                    # Westward
                    P[m-1, m-2] = abs(0.25*pi - alpha + 1.5*pi) / 0.25*pi
            except:
                print i, j, m, n
    return P.tocsr()


def get_upslope_area(R, res=[20.0, 20.0]):

    # Upstream calculation
    def upstream(R, area, i, j):
        """ Recursive function for computing upstream area at (i,j),
        given an array to write to and a flow field.
        """
        def proportion(alpha, position):
            """ Return proportion flowing to center from position given
            flow angle alpha, and where position is
                            0   1   2
                            3   *   4
                            5   6   7
            """
            pi = math.pi

            # Beta is the direction to the centre node
            BETA = [0.75, 1.0, 1.25, 0.5, 1.5, 0.25, 0.0, 1.75]
            beta = BETA[position] * pi

            theta = abs(beta-alpha)
            return theta <= 0.25*pi and (0.25*pi-theta) / (0.25*pi) or 0.0

        if not np.isnan(UA[i,j]):

            return UA[i,j]

        else:

            cumarea = area

            I = [i-1, i-1, i-1,
                 i,        i,
                 i+1, i+1, i+1]
            J = [j-1, j, j+1,
                 j-1,    j+1,
                 j-1, j, j+1]

            try:
                Rij = [R[i-1, j-1],   R[i-1, j],     R[i-1, j+1],
                       R[i, j-1],                     R[i, j+1],
                       R[i+1, j-1],   R[i+1, j],     R[i+1, j+1]]
            except IndexError:
                return 0.0

            iRij = enumerate(Rij)

            # Retains cells where p > 0.0
            contributors = filter(
                    lambda a: proportion(a[1], a[0]) > 0.0, iRij)

            # For each contributor, add the upslope area
            for a in contributors:
                try:
                    cumarea += proportion(a[1], a[0]) * upstream(R, area, I[a[0]], J[a[0]])
                except IndexError:
                    pass

            UA[i,j] = cumarea

            return cumarea


    # Calculate upstream areas
    print "Calculating upstream areas..."
    global UA

    UA = np.nan * np.zeros_like(R)
    area = res[0] * res[1]

    for i in range(1, R.shape[0]):
        for j in range(1, R.shape[1]):
            if np.isnan(UA[i,j]):
                upstream(R, area, i, j)

    return UA
