""" Cython versions of low-level functions for computing vector geometry. """

import numpy as np
cimport numpy as np
from cpython cimport bool
from libc.math cimport sqrt, atan2

def isleft(tuple pt0, tuple pt1, tuple pt2):
    return (pt1[0]-pt0[0])*(pt2[1]-pt0[1]) - (pt1[1]-pt0[1])*(pt2[0]-pt0[0]) > 0.0

def polarangle(tuple pt0, tuple pt1):
    """ Return the polar angle from pt0 to pt1 """
    return atan2(pt1[1]-pt0[1], pt1[0]-pt0[0])

cdef struct Point:
    double x
    double y

cdef struct Vector2:
    double x
    double y

cdef double cdot2(Vector2 u, Vector2 v) nogil:
    return u.x * v.x + u.y * v.y

cdef double cross2(Vector2 u, Vector2 v) nogil:
    return u.x * v.y - u.y * v.x

cdef Point cproj2(Vector2 u, Vector2 v):
    cdef double uv, vv
    uv = cdot2(u, v)
    vv = cdot2(v, v)
    return Point(uv/vv*v.x, uv/vv*v.y)

cdef double dist2(Point pt0, Point pt1) nogil:
    return sqrt((pt0.x-pt1.x)**2 + (pt0.y-pt1.y)**2)

cdef double mind(double a, double b) nogil:
    if a <= b:
        return a
    else:
        return b

cdef double maxd(double a, double b) nogil:
    if a >= b:
        return a
    else:
        return b

cdef double absd(double a) nogil:
    if a < 0:
        return -a
    else:
        return a

def pt_nearest_planar(double x, double y,
                      double endpt0_0, double endpt0_1,
                      double endpt1_0, double endpt1_1):
    """ Determines the point on a segment defined by pairs endpt1
    and endpt2 nearest to a point defined by x, y.
    Returns the a nested tuple of (point, distance)
    """

    cdef double u0, u1, v0, v1
    cdef Point u_on_v, u_int, pt, pt0, pt1
    cdef Vector2 u, v

    pt = Point(x, y)
    pt0 = Point(endpt0_0, endpt0_1)
    pt1 = Point(endpt1_0, endpt1_1)

    u = Vector2(x - pt0.x, y - pt0.y)
    v = Vector2(pt1.x - pt0.x, pt1.y - pt0.y)
    u_on_v = cproj2(u, v)
    u_int = Point(u_on_v.x + pt0.x, u_on_v.y + pt0.y)

    # Determine whether u_int is inside the segment
    # Otherwise return the nearest endpoint
    if absd(pt0.x - pt1.x) < 1e-4:  # Segment is near vertical
        if u_int.y > maxd(pt0.y, pt1.y):
            if pt0.y > pt1.y:
                return ((pt0.x, pt0.y), dist2(pt0, pt))
            else:
                return ((pt1.x, pt1.y), dist2(pt1, pt))

        elif u_int.y < mind(pt0.y, pt1.y):
            if pt0.y < pt1.y:
                return ((pt0.x, pt0.y), dist2(pt0, pt))
            else:
                return ((pt1.x, pt1.y), dist2(pt1, pt))

        else:
            return ((u_int.x, u_int.y), dist2(u_int, pt))

    else:

        if u_int.x > maxd(pt0.x, pt1.x):
            if pt0.x > pt1.x:
                return ((pt0.x, pt0.y), dist2(pt0, pt))
            else:
                return ((pt1.x, pt1.y), dist2(pt1, pt))

        elif u_int.x < mind(pt0.x, pt1.x):
            if pt0.x < pt1.x:
                return ((pt0.x, pt0.y), dist2(pt0, pt))
            else:
                return ((pt1.x, pt1.y), dist2(pt1, pt))

        else:
            return ((u_int.x, u_int.y), dist2(u_int, pt))

# These types could be used for Cython-wrapped Python functions in the future
ctypedef tuple (*fwd_t)(double, double, double, double)
ctypedef tuple (*inv_t)(double, double, double, double)

cdef double _along_distance(object fwd, object inv, double x0, double y0,
        double xp, double yp, double az, double f):
    """ Return the distance between a point distance f along azimuth az away
    from (x0, y0) from (xp, yp). """
    cdef double trialx, trialy, _, d
    (trialx, trialy, _) = fwd(x0, y0, az, f)
    (_, _, d) = inv(trialx, trialy, xp, yp)
    return d

cdef double _along_distance_gradient(object fwd, object inv, double x0, double y0,
        double xp, double yp, double az, double f, double dx):
    """ Return the numerical gradient in terms of f of _along_distance """
    cdef double d1, d2
    d1 = _along_distance(fwd, inv, x0, y0, xp, yp, az, f)
    d2 = _along_distance(fwd, inv, x0, y0, xp, yp, az, f+dx)
    return (d2-d1)/dx

def pt_nearest_proj(object fwd, object inv, tuple pt, double[:] endpt0, double[:] endpt1,
        float tol=0.1, int maxiter=100):
    """ Given geodetic functions *fwd* and *inv*, a Point *pt*, and an arc from
    *endpt1* to *endpt2*, return the point on the arc that is nearest *pt*.

    Scheme employs a bisection minimization method. Iteration continues until a
    tolerance *tol* in meters is reached, or *maxiter* iterations are
    exhausted. If the iteration limit is reached, a ConvergenceError is raised.
    """
    cdef double az, az2, L
    (az, az2, L) = inv(endpt0[0], endpt0[1], endpt1[0], endpt1[1])

    # Detect whether the nearest point is at an endpoint
    cdef double grad0, grad1, d
    grad0 = _along_distance_gradient(fwd, inv, endpt0[0], endpt0[1], pt[0], pt[1], az, 0, 1e-7*L)
    grad1 = _along_distance_gradient(fwd, inv, endpt0[0], endpt0[1], pt[0], pt[1], az, L, 1e-7*L)
    if grad0 > 0:
        d = _along_distance(fwd, inv, endpt0[0], endpt0[1], pt[0], pt[1], az, 0)
        return endpt0, d
    elif grad1 < 0:
        d = _along_distance(fwd, inv, endpt0[0], endpt0[1], pt[0], pt[1], az, L)
        return endpt1, d

    # Bisection iteration
    cdef double x0, x1, xm, xn, yn
    cdef int i = 0
    cdef double dx = tol + 1.0
    cdef double grad
    x0, x1 = 0.0, 1.0
    while dx > tol:
        if i == maxiter:
            raise ConvergenceError("Maximum iterations exhausted in bisection "
                                   "method.")
        xm = 0.5 * (x0 + x1)
        grad = _along_distance_gradient(fwd, inv, endpt0[0], endpt0[1], pt[0], pt[1], az, xm*L, 1e-7*L)
        if grad > 0:
            dx = absd(x1-xm) * L
            x1 = xm
        else:
            dx = absd(x0-xm) * L
            x0 = xm
        i += 1
    (xn, yn, _) = fwd(endpt0[0], endpt0[1], az, xm*L)
    d = _along_distance(fwd, inv, endpt0[0], endpt0[1], pt[0], pt[1], az, xm*L)
    return (xn, yn), d

class ConvergenceError(Exception):
    def __init__(self, message=''):
        self.message = message
    def __str__(self):
        return self.message

# Tree functions
def bbox_intersection_area(tuple bb0, tuple bb1):
    cdef float dx, dy
    dx = maxd(mind(bb0[2], bb1[2]) - maxd(bb0[0], bb1[0]), 0.0)
    dy = maxd(mind(bb0[3], bb1[3]) - maxd(bb0[1], bb1[1]), 0.0)
    return dx*dy

def iswithin(tuple bbox, tuple pt):
    """ Return whether a point is within a bounding box (planar approximation). """
    return (bbox[0] <= pt[0] < bbox[2] and bbox[1] <= pt[1] < bbox[3])

# QuadTree-specific
def hashpt(float xmin, float ymin, float xmax, float ymax, float x, float y):
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
