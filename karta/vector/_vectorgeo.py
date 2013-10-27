""" Low-level functions for computing vector geometry. """

import math
import numpy as np

def intersects(x0, x1, x2, x3, y0, y1, y2, y3):
    """ Determines whether two line segments intersect on a plane. """
    if x1 != x0:
        m0 = float(y1-y0) / float(x1-x0)
    else:
        m0 = 1e37
    if x3 != x2:
        m1 = float(y3-y2) / float(x3-x2)
    else:
        m1 = 1e37
    if m0 == m1:
        return False
    x = float(y1 - y3 + m1*x3 - m0*x1) / float(m1 - m0)
    if (x < max(x0, x1) and x < max(x2, x3) and
        x > min(x0, x1) and x > min(x2, x3)):
        return True
    else:
        return False

def intersections(x0, x1, x2, x3, y0, y1, y2, y3):
    """ Return the point of intersection between two line segments. Returns NaN
    if the lines do not intersect.
    """
    if x1 != x0:
        m0 = float(y1-y0) / float(x1-x0)
    else:
        m0 = 1e37
    if x3 != x2:
        m1 = float(y3-y2) / float(x3-x2)
    else:
        m1 = 1e37
    if m0 == m1:
        return (np.nan, np.nan)
    x = float(m0*x0 - m1*x2 + y2 - y0) / float(m0 - m1)
    if (x < max(x0, x1) and x < max(x2, x3) and
        x > min(x0, x1) and x > min(x2, x3)):
        y = m0 * (x-x0) + y0
        return (x, y)
    else:
        return (np.nan, np.nan)

def distance(pt0, pt1):
    """ Calculate the distance between two points (tuples) """
    d = math.sqrt(sum([abs(a-b)**2 for a, b in zip(pt0, pt1)]))
    return d

def pt_nearest(pt, endpt1, endpt2):
    """ Determines the point on a segment defined by tuples endpt1
    and endpt2 nearest to a point defined by tuple pt.
    Returns the a nested tuple of (point, distance)
    """
    dot2 = lambda v1,v2: float(v1[0]*v2[0] + v1[1]*v2[1])
    proj2 = lambda u,v: [dot2(u,v) / dot2(v,v) * e for e in v]
    dist = lambda u,v: math.sqrt((u[0]-v[0])**2. + (u[1]-v[1])**2.)

    u = (pt[0] - endpt1[0], pt[1] - endpt1[1])
    v = (endpt2[0] - endpt1[0], endpt2[1] - endpt1[1])
    u_on_v = proj2(u,v)
    u_int = (u_on_v[0] + endpt1[0], u_on_v[1] + endpt1[1])
    #dist_u_int = dist(u_int, pt)

    # Determine whether u_int is inside the segment
    # Otherwise return the nearest endpoint
    if endpt1[0] == endpt2[0]:  # Segment is vertical
        if u_int[1] > max(endpt1[1], endpt2[1]):
            return (endpt1[1] > endpt2[1]
                    and (endpt1, dist(endpt1, pt))
                    or (endpt2, dist(endpt2, pt)))
        elif u_int[1] < min(endpt1[1], endpt2[1]):
            return (endpt1[1] < endpt2[1]
                    and (endpt1, dist(endpt1, pt))
                    or (endpt2, dist(endpt2, pt)))
        else:
            return (u_int, dist(u_int, pt))
    else:
        if u_int[0] > max(endpt1[0], endpt2[0]):
            return (endpt1[0] > endpt2[0]
                    and (endpt1, dist(endpt1, pt))
                    or (endpt2, dist(endpt2, pt)))
        elif u_int[0] < min(endpt1[0], endpt2[0]):
            return (endpt1[0] < endpt2[0]
                    and (endpt1, dist(endpt1, pt))
                    or (endpt2, dist(endpt2, pt)))
        else:
            return (u_int, dist(u_int, pt))


