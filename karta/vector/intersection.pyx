import numpy as np
cimport numpy as np
from coordstring cimport CoordString
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

cdef inline double cross_prod2(double u0, double u1, double v0, double v1):
    return u0*v1 - u1*v0

def all_intersections(CoordString a, CoordString b):
    cdef int na = len(a)
    cdef int nb = len(b)
    cdef int i = 0, j = 0
    cdef double x0, x1, x2, x3, y0, y1, y2, y3
    cdef xi, yi
    cdef list intersections = []

    if not a.ring:
        na -= 1
    if not b.ring:
        nb -= 1

    for i in range(na):
        x0 = a.getX(i)
        x1 = a.getX(i+1)
        y0 = a.getY(i)
        y1 = a.getY(i+1)
        for j in range(nb):
            x2 = b.getX(j)
            x3 = b.getX(j+1)
            y2 = b.getY(j)
            y3 = b.getY(j+1)
            xi, yi = intersection(x0, x1, x2, x3, y0, y1, y2, y3)
            if not np.isnan(xi):
                intersections.append((xi, yi))
    return intersections

cpdef intersection(double x0, double x1, double x2, double x3,
                   double y0, double y1, double y2, double y3):
    """ Returns coordinates of intersection point, or (NaN, NaN) if lines don't
    intersect.

    Line1 consists of point pairs 0, 1
    Line2 consists of point pairs 2, 3
    """
    cdef double rxs = cross_prod2(x1-x0, y1-y0, x3-x2, y3-y2)
    if rxs == 0:
        # parallel or collinear
        return np.nan, np.nan

    cdef double t = cross_prod2(x2-x0, y2-y0, x3-x2, y3-y2) / rxs
    cdef double u = cross_prod2(x2-x0, y2-y0, x1-x0, y1-y0) / rxs
    if (0 < t <= 1) and (0 < u <= 1):
        return x0 + t*(x1-x0), y0 + t*(y1-y0)
    else:
        # non-intersecting
        return np.nan, np.nan

def count_crossings(double xp, double yp, CoordString coords):
    """ Count the number of times a vertical ray from (xp, yp) crosses a
    line/polygon defined by *coords*
    """
    cdef int n = len(coords)
    cdef int i, cnt = 0
    if not coords.ring:
        n -= 1
    cdef double x0 = coords.getX(0)
    cdef double y0 = coords.getY(0)
    cdef double x1, y1
    for i in range(1, n):
        x1 = coords.getX(i)
        y1 = coords.getY(i)
        if intersects_cn(xp, yp, x0, x1, y0, y1) == 1:
            cnt += 1
        x0 = x1
        y0 = y1
    return cnt

cdef int intersects_cn(double xp, double yp, double x0, double x1, double y0, double y1):
    """ Test whether a vertical ray emanating up from a point (xp, yp) crosses
    a line segment [(x0, y0), (x1, y1)]. Used to implement a crossing number
    membership test.
    """
    cdef double m, y

    if x0 != x1:
        m = (y1-y0) / (x1-x0)
    else:
        return 0

    y = y0 + m * (xp - x0)

    if y < yp:
        return 0

    cdef int iswithinx = 0
    cdef int iswithiny = 0

    if m > 0.0 and isbetween_incr(y0, y, y1):
        iswithiny = 1
    elif isbetween_incl(y0, y, y1):
        iswithiny = 1
    elif abs(y0-y1) < 1e-15 and abs(y-y0) < 1e-15:
        iswithiny = 1

    if isbetween_incr(x0, xp, x1):
        iswithinx = 1

    return iswithinx*iswithiny

