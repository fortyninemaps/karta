""" Cython versions of low-level functions for computing vector geometry. """

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
        x > dbl_min(x0, y1) and x > dbl_min(x2, x3)):
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
        x > dbl_min(x0, y1) and x > dbl_min(x2, x3)):
        y = m0 * (x-x0) + y0
        return (x, y)
    else:
        return (np.nan, np.nan)
