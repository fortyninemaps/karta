"""
Flow-routing algorithms
"""

import raster
import math
import numpy as np
import pdb

def dinfinity(D, res=(30.0,30.0)):
    """ Infinite flow direction scheme of Tarboton (1997). """
    dx = res[0]
    dy = res[1]

    def get_facet_slope(e0, e1, e2, d1, d2):
        """ Calculate facet slope, where e0 is the central node, e1 is
        a node in a cardinal direction, and e2 is a node in a diagonal.
        """
        s1 = (e0 - e1) / d1
        s2 = (e1 - e2) / d2
        aspect = math.atan(s2 / s1)
        magnitude = math.sqrt(s1*s1 + s2*s2)
        return magnitude, aspect

    def get_flow_dir(A, dx, dy):
        """ For a 3x3 sector A, compute the flow direction, or return
        nan to indicate an unresolved case.
        """
        cells1 = [(0,0), (0,1), (0,2),
                  (1,0),        (1,2),
                  (2,0), (2,1), (2,2)]
        cells2 = [(0,1), (0,2), (1,2),
                  (0,0),        (2,2),
                  (1,0), (2,0), (2,1)]
        d1 = [dx, dx, dy,
              dy,     dy,
              dy, dx, dx]
        d2 = math.sqrt(dx*dx + dy*dy)

        # Compute all slopes
        fs = map(get_facet_slope, [A[1,1] for c in cells1],
                                  [A[c] for c in cells1],
                                  [A[c] for c in cells2],
                                  d1, [d2 for i in d1])

        # Return the direction of the most negative, or return nan
        pdb.set_trace()
        fs.sort()
        fs.reverse()
        return fs[0][0] < 0 and fs[0][1] or np.nan


    # For each cell in D[1:-1, 1:-1], find the flow direction
    #   First, associate a 3x3 matrix with each cell, composed of neighbouring cells
    #   Next, apply a vectorized routine to do the calculation

    D3 = np.array([
                  [D[:-2,:-2], D[1:-1,:-2], D[2:,:-2]],
                  [D[:-2,1:-1], D[1:-1,1:-1], D[2:,1:-1]],
                  [D[:-2,2:], D[1:-1,2:], D[2:,2:]],
                  ]).transpose([2,3,0,1])

    #get_flow_dir_v = np.vectorize(get_flow_dir)

    flow_dir = np.zeros_like(D[1:-1,1:-1])
    for i in range(D3.shape[0]):
        for j in range(D3.shape[1]):
            flow_dir[i,j] = get_flow_dir(D3[i,j,:,:], dx, dy)

    return flow_dir         # temporary

    #flow_dir = get_flow_dir_v(D3, dx*np.ones_like(D[1:-1,1:-1]),
                                  #dy*np.ones_like(D[1:-1,1:-1]))


    # Upstream calculation
    def upstream(fd, area, i, j):
        """ Recursive function for computing upstream area at (i,j),
        given an array to write to and a flow field.
        """

        def proportion(alpha, position):
            """ Return proportion flowing to center from position given
            angle alpha, and where position is
                            0   1   2
                            3   *   4
                            5   6   7
            """
            pi = math.pi

            if position == 0:
                pass
            elif position == 1:
                alpha -= 0.25*pi
            elif position == 2:
                alpha -= 0.5*pi
            elif position == 3:
                alpha += 0.25*pi
            elif position == 4:
                alpha += 0.75*pi
            elif position == 5:
                alpha -= 0.5*pi
            elif position == 6:
                alpha -= 0.75*pi
            elif position == 7:
                alpha += pi

            if alpha > 0.25*pi:
                alpha = 0.5*pi - alpha
            if alpha < 0:
                alpha = 0

            return alpha / 0.25*pi

        if UA[i,j] is not np.nan:
            return UA[i,j]

        else:
            cumarea = area

            fdij = [fd[i-1, j-1],   fd[i-1, j],     fd[i-1, j+1],
                    fd[i, j-1],                     fd[i, j+1],
                    fd[i+1, j-1],   fd[i+1, j],     fd[i+1, j+1]]

            # Test whether the surrounding cells flow into (i,j)
            for position, angle in enumerate(fdij):
                cumarea += proportion(angle, position) * upstream(fd, k, l)

            return cumarea

    # Calculate upstream areas
    global UA

    UA = np.nan * np.zeros_like(D)
    area = res[0] * res[1]

    for i in range(1, D.shape[0]-1):
        for j in range(1, D.shape[1]-1):

            if UA[i,j] is not np.nan:

                upstream(flow_dir, area, i, j)

    return UA






def d8_mult(D):
    """ Multiple flow direction D8 algorithm.

    Described by Quinn et al. (1991).
    """
    F = np.zeros_like(D)
    dxy = math.sqrt(dx8dx+dy*dy)

    # Calculate gradients
    dzdx = (D[:,1:] - D[:,:-1]) / dx
    dzdy = (D[1:,:] - D[:-1,:]) / dy
    dzd1 = (D[1:,1:] - D[:-1,:-1]) / dxy
    dzd2 = (D[1:,:-1] - D[:-1,1:]) / dxy

    # Define a recursive function to calculate upstream area
    def upstream():
        pass


    # Vectorize
    upstreamv = np.vectorize(upstream)

    # Apply to D
    F = upstreamv(D)



    return F




