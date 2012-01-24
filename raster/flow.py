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
import pdb, time

def dinfinity(D, res=(30.0,30.0)):
    """ Infinite flow direction scheme of Tarboton (1997). """
    t0 = time.time()
    dx = res[0]
    dy = res[1]

    def get_facet_slope_old(e0, e1, e2, d1, d2, ac, af):
        #""" Calculate facet slope, where e0 is the central node, e1 is
        #a node in a cardinal direction, and e2 is a node in a diagonal.
        #"""
        ##print e0,e1,e2, "\t", d1,d2, "\t", ac,af
        #pi = math.pi
        #s1 = (e0 - e1) / d1
        #s2 = (e1 - e2) / d2
        #try:
            #aspect = math.atan(s2 / s1)
        #except ZeroDivisionError:
            #aspect = 0.25 * pi
        #magnitude = math.sqrt(s1*s1 + s2*s2)

        #if aspect > 0.25*pi:
            #aspect = 0.25*pi
            #magnitude = abs((e0 - e2) / math.sqrt(d1*d1 + d2*d2))
        #elif aspect < -0.25*pi:
            #aspect = -0.25*pi
            #magnitude = abs((e0 - e2) / math.sqrt(d1*d1 + d2*d2))

        ## Rotate / flip facet accordng to ac, af
        #aspect_corr = af * aspect + ac * 0.5*pi
        ##aspect_corr = aspect

        ##print "\t", magnitude, aspect_corr
        #return magnitude, aspect_corr
        return

    def get_facet_slope(e0, e1, e2, d1, d2, ac, af):
        """ Calculate facet slope, where e0 is the central node, e1 is
        a node in a cardinal direction, and e2 is a node in a diagonal.
        """
        pi = math.pi

        u1 = (0.0, d1, e1-e0)
        u2 = (d2, 0.0, e2-e1)
        n = -np.cross(u1, u2)
        nhat = n / np.linalg.norm(n)
        nhat[nhat==0.0] = 0.0
        aspect = math.atan(nhat[0] / nhat[1])
        slope = math.acos(nhat[2])

        if nhat[1] < 0:
            aspect += pi

        # For aspects that point inward, flip the sign of slope
        if aspect < 0:
            aspect += 2.0*pi
        if aspect > 0.5*pi and aspect <= 1.5*pi:
            slope = -slope

        # Rotate / flip facet accordng to ac, af
        aspect_corr = af * aspect + ac * 0.5*pi

        # Make sure aspect is in the range [0, 2pi)
        while aspect_corr >= 2*pi:
            aspect_corr -= 2*pi
        while aspect_corr < 0.0:
            aspect_corr += 2*pi

        return slope, aspect_corr

    def get_flow_dir(A, dx, dy):
        """ For a 3x3 sector A, compute the direction of maximum outward
        slope (assumed to be the flow direction). Returns nan when
        indeterminate.
        """
        e1 = [(0,1), (0,1), (1,2),
              (1,0),        (1,2),
              (1,0), (2,1), (2,1)]
        e2 = [(0,0), (0,2), (0,2),
              (0,0),        (2,2),
              (2,0), (2,0), (2,2)]
        d1 = [dy, dy, dx,
              dx,     dx,
              dx, dy, dy]
        d2 = [dx, dx, dy,
              dy,     dy,
              dy, dx, dx]
        ac = [0, 0, 1,              # Translation table
              -1,   1,
              -1, -2, 2]
        af = [-1, 1, -1,            # Flipping table
              1,      1,
              -1, 1, -1]

        # Compute all slopes
        fs = map(get_facet_slope, [A[1,1] for _e in e1],
                                  [A[_e] for _e in e1],
                                  [A[_e] for _e in e2],
                                  d1, d2, ac, af)

        # Return the direction of the steepest inward slope, or return 0.0
        fs.sort()
        return fs[0][0] < 0 and fs[0][1] or np.nan


    # For each cell in D[1:-1, 1:-1], find the flow direction
    print "Calculating flow directions..."

    # Associate a 3x3 matrix with each cell, composed of neighbouring cells
    D3 = np.array([
                  [D[:-2,:-2], D[1:-1,:-2], D[2:,:-2]],
                  [D[:-2,1:-1], D[1:-1,1:-1], D[2:,1:-1]],
                  [D[:-2,2:], D[1:-1,2:], D[2:,2:]],
                  ]).transpose([2,3,1,0])


    flow_field = np.zeros_like(D[1:-1,1:-1])
    for i in range(D3.shape[0]):
        for j in range(D3.shape[1]):
            flow_field[i,j] = get_flow_dir(D3[i,j,:,:], dx, dy)

    #get_flow_dir_v = np.vectorize(get_flow_dir)
    #flow_field = get_flow_dir_v(D3, dx*np.ones_like(D[1:-1,1:-1]),
                                  #dy*np.ones_like(D[1:-1,1:-1]))

    print "\tdone"
    #return flow_field         # temporary

    # Upstream calculation
    def upstream(ff, area, i, j):
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
                ffij = [ff[i-1, j-1],   ff[i-1, j],     ff[i-1, j+1],
                        ff[i, j-1],                     ff[i, j+1],
                        ff[i+1, j-1],   ff[i+1, j],     ff[i+1, j+1]]
            except IndexError:
                return 0.0

            iffij = enumerate(ffij)

            # Retains cells where p > 0.0
            contributors = filter(
                    lambda a: proportion(a[1], a[0]) > 0.0, iffij)

            # For each contributor, add the upslope area
            for a in contributors:
                try:
                    cumarea += proportion(a[1], a[0]) * upstream(ff, area, I[a[0]], J[a[0]])
                except IndexError:
                    pass

            UA[i,j] = cumarea

            return cumarea


    # Calculate upstream areas
    print "Calculating upstream areas..."
    global UA

    UA = np.nan * np.zeros_like(D)
    area = res[0] * res[1]

    for i in range(1, D.shape[0]):
        for j in range(1, D.shape[1]):
            if np.isnan(UA[i,j]):
                upstream(flow_field, area, i, j)

    print "\tdone"

    print "time:", time.time() - t0

    return raster.pad(UA)


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

