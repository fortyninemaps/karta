""" Cython versions of low-level functions for computing vector geometry. """

import math
import numpy as np
cimport numpy as np

cdef inline double dbl_max(double a, double b): return a if a >= b else b
cdef inline double dbl_min(double a, double b): return a if a <= b else b

cpdef int intersects(double x0, double x1, double x2, double x3,
                     double y0, double y1, double y2, double y3):
    """ Determines whether two line segments intersect on a plane. The line
    segments are specified as ((x0,y0),(x1,y1)), etc."""
    cdef double m0, m1
    cdef double x

    if x1 != x0:
        m0 = (y1-y0) / (x1-x0)
    else:
        m0 = 1e37
    if x3 != x2:
        m1 = (y3-y2) / (x3-x2)
    else:
        m1 = 1e37
    if m0 == m1:
        return False
    x = (y1 - y3 + m1*x3 - m0*x1) / (m1 - m0)
    if (x < dbl_max(x0, x1) and x < dbl_max(x2, x3) and
        x > dbl_min(x0, x1) and x > dbl_min(x2, x3)):
        return True
    else:
        return False

def intersections(double x0, double x1, double x2, double x3,
                 double y0, double y1, double y2, double y3):
    """ Return the point of intersection between two line segments. Returns NaN
    if the lines do not intersect.
    """
    cdef double m0, m1
    cdef double x, y

    if x1!=x0:
        m0 = (y1-y0) / (x1-x0)
    else:
        m0 = 1e16
    if x3!=x2:
        m1 = (y3-y2) / (x3-x2)
    else:
        m1 = 1e16
    if m0 == m1:
        return (np.nan, np.nan)
    x = (m0*x0 - m1*x2 + y2 - y0) / (m0 - m1)
    if (x < dbl_max(x0, x1) and x < dbl_max(x2, x3) and
        x > dbl_min(x0, x1) and x > dbl_min(x2, x3)):
        y = m0 * (x-x0) + y0
        return (x, y)
    else:
        return (np.nan, np.nan)

cdef double cdot2(tuple v1, tuple v2):
    return v1[0] * v2[0] + v1[1] * v2[1]

cdef tuple cproj2(tuple u, tuple v):
    cdef double uv, vv
    uv = cdot2(u,v)
    vv = cdot2(v,v)
    return (uv / vv * v[0], uv / vv * v[1])

def distance(tuple pt0, tuple pt1):
    """ Calculate the distance between two points (tuples) """
    cdef float dsq
    dsq = 0.0
    for i in range(min(len(pt0), len(pt1))):
        dsq += abs(pt0[i] - pt1[i])**2
    return math.sqrt(dsq)

def pt_nearest(tuple pt, tuple endpt1, tuple endpt2):
    """ Determines the point on a segment defined by tuples endpt1
    and endpt2 nearest to a point defined by tuple pt.
    Returns the a nested tuple of (point, distance)
    """
    
    cdef tuple u, v, u_int
    cdef tuple u_on_v
    
    dist = lambda u,v: math.sqrt((u[0]-v[0])**2. + (u[1]-v[1])**2.)

    u = (pt[0] - endpt1[0], pt[1] - endpt1[1])
    v = (endpt2[0] - endpt1[0], endpt2[1] - endpt1[1])
    u_on_v = cproj2(u,v)
    u_int = (u_on_v[0] + endpt1[0], u_on_v[1] + endpt1[1])

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


