import numpy as np
cimport numpy as np
from cpython cimport bool

cdef inline double dbl_max(double a, double b): return a if a >= b else b
cdef inline double dbl_min(double a, double b): return a if a <= b else b

cdef bool isbetween_inc(double a, double b, double c):
    return dbl_min(a, c) <= b <= dbl_max(a, c)

cdef bool isbetween_incl(double a, double b, double c):
    return dbl_min(a, c) <= b < dbl_max(a, c)

cdef bool isbetween_incr(double a, double b, double c):
    return dbl_min(a, c) < b <= dbl_max(a, c)

ctypedef bool (*isbetween_t)(double, double, double)

def intersection(double x0, double x1, double x2, double x3,
                 double y0, double y1, double y2, double y3):
    """ Return the point of intersection between two line segments on a plane.
    Returns (NaN, NaN) if the lines do not intersect. """
    cdef double m0, m1
    cdef double x, y

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

    cdef bool iswithinx
    cdef bool iswithiny
    iswithinx = False
    iswithiny = False

    if abs(x1 - x0) >= 1e-15:
        if abs(x3 - x2) >= 1e-15:
            if isbetween_inc(x0, x, x1) and isbetween_inc(x2, x, x3):
                iswithinx = True
        else:
            if abs(x - x2) < 1e-15 and isbetween_inc(x0, x, x1):
                iswithinx = True
    else:
        if abs(x - x0) < 1e-15 and isbetween_inc(x2, x, x3):
            iswithinx = True

    if iswithinx:
        if abs(x-x0) >= 1e-15:
            y = m0 * (x-x0) + y0
        else:
            y = m1 * (x-x2) + y2
        if abs(y1 - y0) >= 1e-15:
            if abs(y3 - y2) >= 1e-15:
                if isbetween_inc(y0, y, y1) and isbetween_inc(y2, y, y3):
                    iswithiny = True
            elif abs(y - y2) < 1e-15 and isbetween_inc(y0, y, y1):
                iswithiny = True
        elif abs(y - y0) < 1e-15 and isbetween_inc(y2, y, y3):
            iswithiny = True

    if iswithinx and iswithiny:
        return (x, y)
    else:
        return (np.nan, np.nan)

def intersects_cn(double xp, double yp, double x0, double x1, double y0, double y1):
    """ Test whether a vertical ray emanating up from a point (xp, yp) crosses
    a line segment [(x0, y0), (x1, y1)]. Used to implement a crossing number
    membership test.
    """
    cdef double m, y

    if x0 != x1:
        m = (y1-y0) / (x1-x0)
    else:
        return False

    y = y0 + m * (xp - x0)

    if y < yp:
        return False

    cdef bool iswithinx
    cdef bool iswithiny
    iswithinx = False
    iswithiny = False

    if m > 0.0 and isbetween_incr(y0, y, y1):
        iswithiny = True
    elif isbetween_incl(y0, y, y1):
        iswithiny = True
    elif abs(y0-y1) < 1e-15 and abs(y-y0) < 1e-15:
        iswithiny = True

    if isbetween_incr(x0, xp, x1):
        iswithinx = True

    return iswithinx and iswithiny

