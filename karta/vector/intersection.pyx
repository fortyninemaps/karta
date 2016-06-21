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
    for i in range(na-1):
        x0 = a.coords[i*a.rank]
        x1 = a.coords[(i+1)*a.rank]
        y0 = a.coords[i*a.rank+1]
        y1 = a.coords[(i+1)*a.rank+1]
        for j in range(nb-1):
            x2 = b.coords[j*b.rank]
            x3 = b.coords[(j+1)*b.rank]
            y2 = b.coords[j*b.rank+1]
            y3 = b.coords[(j+1)*b.rank+1]
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
    if (1e-14 <= t <= 1) and (1e-14 <= u <= 1):
        return x0 + t*(x1-x0), y0 + t*(y1-y0)
    else:
        # non-intersecting
        return np.nan, np.nan

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

