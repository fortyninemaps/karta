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

    #def get_facet_slope(e0, e1, e2, d1, d2, ac, af):
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

        # For aspects that point backward, flip the sign of slope
        #if aspect > 0.5*pi and aspect <= 1.5*pi:
            #slope = -slope

        # Rotate / flip facet accordng to ac, af
        aspect_corr = af * aspect + ac * 0.5*pi

        # Make sure aspect is in the range [0, 2pi)
        while aspect_corr >= 2*pi:
            aspect_corr -= 2*pi

        return slope, aspect_corr

    def get_flow_dir(A, dx, dy):
        """ For a 3x3 sector A, compute the flow direction, or return
        nan to indicate an unresolved case.
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

        FS = np.array([[fs[0][0], fs[1][0], fs[2][0]],
                       [fs[3][0], 0.0,   fs[4][0]],
                       [fs[5][0], fs[6][0], fs[7][0]]])
        pi = math.pi
        AS = np.array([[fs[0][1]/pi, fs[1][1]/pi, fs[2][1]/pi],
                       [fs[3][1]/pi, 0.0,         fs[4][1]/pi],
                       [fs[5][1]/pi, fs[6][1]/pi, fs[7][1]/pi]])

        # Return the direction of the steepest slope, or return 0.0
        fs.sort()
        fs.reverse()
        return fs[-1][0] > 0 and fs[-1][1] or 0.0


    # For each cell in D[1:-1, 1:-1], find the flow direction
    #   First, associate a 3x3 matrix with each cell, composed of neighbouring cells
    #   Next, apply a vectorized routine to do the calculation

    D3 = np.array([
                  [D[:-2,:-2], D[1:-1,:-2], D[2:,:-2]],
                  [D[:-2,1:-1], D[1:-1,1:-1], D[2:,1:-1]],
                  [D[:-2,2:], D[1:-1,2:], D[2:,2:]],
                  ]).transpose([2,3,1,0])

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




