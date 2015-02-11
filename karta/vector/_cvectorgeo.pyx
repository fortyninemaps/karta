""" Cython versions of low-level functions for computing vector geometry. """

import numpy as np
cimport numpy as np
from cpython cimport bool
from libc.math cimport atan, sqrt

cpdef double cdot2(tuple v1, tuple v2):
    return v1[0] * v2[0] + v1[1] * v2[1]

cpdef double cross2(tuple pt1, tuple pt2):
    return pt1[0]*pt2[1] - pt1[1]*pt2[0]

def isleft(tuple pt0, tuple pt1, tuple pt2):
    return (pt1[0]-pt0[0])*(pt2[1]-pt0[1]) - (pt1[1]-pt0[1])*(pt2[0]-pt0[0]) > 0.0

def polarangle(tuple pt0, tuple pt1):
    """ Return the polar angle from pt0 to pt1, where we assume pt0.y
    <= pt1.y """
    cdef double pi
    cdef double dx
    pi = 3.141592653589793
    dx = pt1[0] - pt0[0]
    if dx > 0:
        return atan((pt1[1] - pt0[1]) / dx)
    elif dx < 0:
        return atan((pt1[1] - pt0[1]) / dx) + pi
    else:
        return 0.5*pi

cdef inline double dbl_max(double a, double b): return a if a >= b else b
cdef inline double dbl_min(double a, double b): return a if a <= b else b

cdef bool isbetween_inc(double a, double b, double c):
    return b <= dbl_max(a, c) and b >= dbl_min(a, c)

cdef bool isbetween_incl(double a, double b, double c):
    return b < dbl_max(a, c) and b >= dbl_min(a, c)

cdef bool isbetween_incr(double a, double b, double c):
    return b <= dbl_max(a, c) and b > dbl_min(a, c)

ctypedef bool (*isbetween_t)(double, double, double)

def intersection(double x0, double x1, double x2, double x3,
                 double y0, double y1, double y2, double y3):
    """ Return the point of intersection between two line segments. Returns NaN
    if the lines do not intersect.
    """
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
    """ Test whether a horizontal ray emanating left from a point (xp, yp)
    crosses a line segment. Used to implement a crossing number membership
    test.
    """
    cdef double m0, m1
    cdef double x, y

    if x1 != x0:
        m = float(y1-y0) / float(x1-x0)
    else:
        m = 1e37
    if m == 0.0:
        return False

    x = float(yp-y0) / m  + x0

    if x < xp:
        return False
    
    cdef bool iswithinx
    cdef bool iswithiny
    iswithinx = False
    iswithiny = False

    cdef isbetween_t isbetween
    if m > 0:
        isbetween = isbetween_incl
    else:
        isbetween = isbetween_incr

    if abs(x1 - x0) >= 1e-15 and isbetween(x0, x, x1):
        iswithinx = True
    elif abs(x - x0) < 1e-15:
        iswithinx = True

    if iswithinx:
        if abs(y1 - y0) >= 1e-15 and isbetween(y0, yp, y1):
            iswithiny = True
        elif abs(yp - y0) < 1e-15:
            iswithiny = True
    return iswithinx and iswithiny

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
    return sqrt(dsq)

def pt_nearest_planar(tuple pt, tuple endpt1, tuple endpt2):
    """ Determines the point on a segment defined by tuples endpt1
    and endpt2 nearest to a point defined by tuple pt.
    Returns the a nested tuple of (point, distance)
    """
    
    cdef tuple u, v, u_int
    cdef tuple u_on_v
    
    dist = lambda u,v: sqrt((u[0]-v[0])**2. + (u[1]-v[1])**2.)

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

def pt_nearest_proj(fwd, inv, pt, endpt0, endpt1, float tol=1.0, int maxiter=50):
    """ Given geodetic functions *fwd* and *inv*, a Point *pt*, and an arc from
    *endpt1* to *endpt2*, return the point on the arc that is nearest *pt*.

    Scheme employs a bisection minimization method. Iteration continues until a
    tolerance *tol* in meters is reached, or *maxiter* iterations are
    exhausted. If the iteration limit is reached, a ConvergenceError is raised.
    """
    cdef double az, az2, L
    (az, az2, L) = inv(endpt0[0], endpt0[1], endpt1[0], endpt1[1])

    def distance(double x):
        cdef double trialx, trialy, _, d
        (trialx, trialy, _) = fwd(endpt0[0], endpt0[1], az, x*L)
        (_, _, d) = inv(trialx, trialy, pt[0], pt[1])
        return d
    
    def ddx(double x):
        cdef double dx, d1, d2
        dx = 1e-8
        d1 = distance(x)
        d2 = distance(x+dx)
        return (d2-d1) / dx
    
    # Detect whether the nearest point is at an endpoint
    cdef double dx0, dx1
    dx0, dx1 = ddx(0), ddx(1)
    if dx0 > 0:
        return endpt0, distance(0)
    elif dx1 < 0:
        return endpt1, distance(1)
    
    # Bisection iteration
    cdef float x0, x1, xm, xn, yn
    cdef int i = 0
    cdef float dx = tol + 1.0
    x0, x1 = 0.0, 1.0
    while dx > tol:
        if i == maxiter:
            raise ConvergenceError("Maximum iterations exhausted in bisection "
                                   "method.")
        xm = 0.5 * (x0 + x1)
        if ddx(xm) > 0:
            dx = abs(x1-xm) * L
            x1 = xm
        else:
            dx = abs(x0-xm) * L
            x0 = xm
        i += 1
    (xn, yn, _) = fwd(endpt0[0], endpt0[1], az, xm*L)
    return (xn, yn), distance(xm)

class ConvergenceError(Exception):
    def __init__(self, message=''):
        self.message = message
    def __str__(self):
        return self.message


# QuadTree functions

def iswithin(tuple bbox, tuple pt):
    """ Return whether a point is within a bounding box (planar approximation). """
    cdef float xmn = bbox[0]
    cdef float xmx = bbox[1]
    cdef float ymn = bbox[2]
    cdef float ymx = bbox[3]
    cdef float x = pt[0]
    cdef float y = pt[1]
    if (xmn <= x < xmx and ymn <= y < ymx):
        return True
    else:
        return False

def hashpt(float xmin, float xmax, float ymin, float ymax, float x, float y):
    """ Returns a generator that returns successive quadrants [0-3] that
    constitute a geohash for *pt* in a global *bbox*. """
    cdef float xm, ym
    cdef int geohash
    while True:
        xm = 0.5 * (xmin + xmax)
        ym = 0.5 * (ymin + ymax)
        if x < xm:
            if y < ym:
                geohash = 0
                xmax = xm
                ymax = ym
            elif y >= ym:
                geohash = 2
                xmax = xm
                ymin = ym
            else:
                raise HashError
        elif x >= xm:
            if y < ym:
                geohash = 1
                xmin = xm
                ymax = ym
            elif y >= ym:
                geohash = 3
                xmin = xm
                ymin = ym
            else:
                raise HashError
        else:
            raise HashError
        yield geohash
        
class HashError(Exception):
    pass

