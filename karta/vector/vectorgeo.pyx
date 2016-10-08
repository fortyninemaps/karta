""" Cython versions of low-level functions for computing vector geometry. """

from libc.math cimport (NAN, M_PI, sqrt,
                        sin, cos, tan, asin, acos, atan, atan2,
                        fmin, fmax, fabs)
from cpython cimport bool
cimport cython
from coordstring cimport CoordString

cdef struct Vector2:
    double x
    double y

cdef struct Vector3:
    double x
    double y
    double z

cdef inline double dot2(Vector2 u, Vector2 v) nogil:
    return u.x*v.x + u.y*v.y

cdef inline double dot3(Vector3 u, Vector3 v) nogil:
    return u.x*v.x + u.y*v.y + u.z*v.z

cdef inline double cross2(Vector2 u, Vector2 v) nogil:
    return u.x*v.y - u.y*v.x

cdef inline Vector3 cross3(Vector3 u, Vector3 v):
    return Vector3(u.y*v.z - u.z*v.y, u.x*v.z - u.z*v.x, u.x*v.y - u.y*v.x)

cdef Vector2 proj2(Vector2 u, Vector2 v):
    cdef double uv_vv
    uv_vv = dot2(u, v) / dot2(v, v)
    return Vector2(uv_vv*v.x, uv_vv*v.y)

cdef inline double dist2(Vector2 pt0, Vector2 pt1) nogil:
    return sqrt((pt0.x-pt1.x)**2 + (pt0.y-pt1.y)**2)

@cython.cdivision(True)
cdef double dist_sph(Vector2 pt0, Vector2 pt1) nogil:
    """ Returns distance on a sphere of radius 1 """
    cdef double dx = fabs(pt0.x - pt1.x) * M_PI / 180.0
    cdef double dy = fabs(pt0.y - pt1.y) * M_PI / 180.0
    cdef double y0 = pt0.y * M_PI / 180.0
    cdef double y1 = pt1.y * M_PI / 180.0
    cdef double d_ = 0.0

    if (dx > 0.01) or (dy > 0.01):
        # use spherical law of cosines
        d_ = acos(sin(y0)*sin(y1) + cos(y0)*cos(y1)*cos(dx))

    else:
        # use haversine
        d_ = (2 * asin(sqrt(sin(0.5*dy)**2 + cos(y0)*cos(y1)*sin(0.5*dx)**2)))
    return d_

cdef int fsign(double a):
    """ Return the sign of *a* """
    if a == 0.0:
        return 1
    else:
        return int(a/fabs(a))

cdef int bndlat_sph(Vector2 pt0, Vector2 pt1, double *ymin, double *ymax):
    """ Return the bounding latitudes of a great circle segment on a sphere.
    Returns 0 on success and 1 if the segment is degenerate. """
    cdef int s0 = fsign(pt0.y)
    cdef int s1 = fsign(pt1.y)
    cdef double faz, baz

    cdef double dlam = (pt1.x-pt0.x) * M_PI / 180.0
    cdef double phi0 = pt0.y * M_PI / 180.0
    cdef double phi1 = pt1.y * M_PI / 180.0

    if dlam != 0.0:
        faz = atan2(sin(dlam)*cos(phi1),
                    cos(phi0)*sin(phi1) - sin(phi0)*cos(phi1)*cos(dlam))
        baz = atan2(sin(-dlam)*cos(phi0),
                    cos(phi1)*sin(phi0) - sin(phi1)*cos(phi0)*cos(-dlam))
    elif phi0 == phi1:
        return 1

    if dlam == 0.0 or (s0 != s1):
        ymin[0] = fmin(pt0.y, pt1.y)
        ymax[0] = fmax(pt0.y, pt1.y)
    else:
        if s0 == 1: # northern
            ymin[0] = fmin(phi0, phi1) * 180.0 / M_PI

            if ((fabs((faz + M_PI) % (2*M_PI) - M_PI) < 0.5*M_PI) and
                (fabs((baz + M_PI) % (2*M_PI) - M_PI) < 0.5*M_PI)):
                # both azimuths "point up"
                ymax[0] = acos(fabs(sin(faz)*cos(phi0))) * 180.0 / M_PI
            else:
                ymax[0] = fmax(phi0, phi1) * 180.0 / M_PI

        else:
            ymax[0] = fmax(phi0, phi1) * 180.0 / M_PI

            if ((fabs((faz + M_PI) % (2*M_PI) - M_PI) > 0.5*M_PI) and
                (fabs((baz + M_PI) % (2*M_PI) - M_PI) > 0.5*M_PI)):
                # both azimuths "point down"
                ymin[0] = -acos(fabs(sin(faz)*cos(phi0))) * 180.0 / M_PI
            else:
                ymin[0] = fmin(phi0, phi1) * 180.0 / M_PI
    return 0

def bbox(CoordString cs):
    cdef double xmin, xmax, ymin, ymax, x, y, x0, y0
    cdef int i = 1, n = len(cs)
    if n == 0:
        return (NAN, NAN, NAN, NAN)
    x0 = cs.getX(0)
    y0 = cs.getY(0)
    xmin = x0
    xmax = x0
    ymin = y0
    ymax = y0
    while i != n:
        x = cs.getX(i)
        y = cs.getY(i)
        xmin = fmin(xmin, x)
        xmax = fmax(xmax, x)
        ymin = fmin(ymin, y)
        ymax = fmax(ymax, y)
        x0 = x
        y0 = y
        i += 1
    return (xmin, ymin, xmax, ymax)

@cython.cdivision(True)
cdef inline Vector3 sph2cart(Vector2 a):
    """ Convert a (lambda, phi) coordinate on a sphere with an origin at
    (0, 0, 0) to an (x, y, z) coordinate. """
    cdef double theta = 90.0 - a.y
    return Vector3(sin(M_PI*theta/180.0)*cos(M_PI*a.x/180.0),
                   sin(M_PI*theta/180.0)*sin(M_PI*a.x/180.0),
                   cos(M_PI*theta/180.0))

@cython.cdivision(True)
cdef inline Vector2 cart2sph(Vector3 a):
    """ Convert an (x, y, z) coordinate to a (lambda, phi) coordinate on a
    sphere with an origin at (0, 0, 0). """
    cdef double lon = 0.0
    cdef double lat = 0.0
    if fabs(a.x) > 1e-8:
        lon = atan2(a.y, a.x)
    else:
        lon = asin(a.y / sqrt(a.x*a.x + a.y*a.y))
    if fabs(a.z) > 1e-8:
        lat = 0.5*M_PI - atan(sqrt(a.x*a.x + a.y*a.y)/a.z)
    else:
        lat = 0.5*M_PI - acos(a.z / sqrt(a.x*a.x + a.y*a.y + a.z*a.z))
    return Vector2(lon*180.0/M_PI, lat*180.0/M_PI)

cdef Vector3 eulerpole_cart(Vector3 a, Vector3 b):
    return cross3(a, b)

cdef Vector3 eulerpole(Vector2 a, Vector2 b):
    cdef Vector3 ac = sph2cart(a)
    cdef Vector3 bc = sph2cart(b)
    cdef ep = cross3(ac, bc)
    return ep

cdef double azimuth(Vector2 pt1, Vector2 pt2):
    """ Return azimuth in radians between two points on a plane. """
    return atan2(pt2.x-pt1.x, pt2.y-pt1.y)

@cython.cdivision(True)
cdef double azimuth_sph(Vector2 pt1, Vector2 pt2):
    """ Return azimuth in radians between two points on a sphere. Inputs taken
    as degrees of (lon, lat). """
    cdef double dlon = (pt2.x - pt1.x) * M_PI / 180.0
    cdef double y1 = pt1.y * M_PI / 180.0
    cdef double y2 = pt2.y * M_PI / 180.0
    return atan2(sin(dlon), cos(y1) * tan(y2) - sin(y1) * cos(dlon))

def length(CoordString cs):
    """ Compute planar length of CoordString """
    cdef int n = len(cs)
    if cs.ring:
        n += 1
    cdef int i = 0
    cdef double d = 0.0
    cdef double x0, y0, x1, y1
    x0 = cs.getX(i)
    y0 = cs.getY(i)
    while i != (n-1):
        x1 = cs.getX(i+1)
        y1 = cs.getY(i+1)
        d += sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))
        x0 = x1
        y0 = y1
        i += 1
    return d

def pt_nearest_planar(double x, double y,
                      double endpt0_0, double endpt0_1,
                      double endpt1_0, double endpt1_1):
    """ Determines the point on a segment defined by pairs endpt1
    and endpt2 nearest to a point defined by x, y.
    Returns the a nested tuple of (point, distance)
    """

    cdef double u0, u1, v0, v1
    cdef Vector2 u_on_v, u_int, pt, pt0, pt1
    cdef Vector2 u, v

    pt = Vector2(x, y)
    pt0 = Vector2(endpt0_0, endpt0_1)
    pt1 = Vector2(endpt1_0, endpt1_1)

    u = Vector2(x - pt0.x, y - pt0.y)
    v = Vector2(pt1.x - pt0.x, pt1.y - pt0.y)
    u_on_v = proj2(u, v)
    u_int = Vector2(u_on_v.x + pt0.x, u_on_v.y + pt0.y)

    # Determine whether u_int is inside the segment
    # Otherwise return the nearest endpoint
    if fabs(pt0.x - pt1.x) < 1e-4:  # Segment is near vertical
        if u_int.y > fmax(pt0.y, pt1.y):
            if pt0.y > pt1.y:
                return ((pt0.x, pt0.y), dist2(pt0, pt))
            else:
                return ((pt1.x, pt1.y), dist2(pt1, pt))

        elif u_int.y < fmin(pt0.y, pt1.y):
            if pt0.y < pt1.y:
                return ((pt0.x, pt0.y), dist2(pt0, pt))
            else:
                return ((pt1.x, pt1.y), dist2(pt1, pt))

        else:
            return ((u_int.x, u_int.y), dist2(u_int, pt))

    else:

        if u_int.x > fmax(pt0.x, pt1.x):
            if pt0.x > pt1.x:
                return ((pt0.x, pt0.y), dist2(pt0, pt))
            else:
                return ((pt1.x, pt1.y), dist2(pt1, pt))

        elif u_int.x < fmin(pt0.x, pt1.x):
            if pt0.x < pt1.x:
                return ((pt0.x, pt0.y), dist2(pt0, pt))
            else:
                return ((pt1.x, pt1.y), dist2(pt1, pt))

        else:
            return ((u_int.x, u_int.y), dist2(u_int, pt))

# Currently unused, but could be employed to optimize pyproj function calls
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
    """ Given geodetic functions *fwd* and *inv*, a *pt*, and an arc from
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
    cdef double x0 = 0.0, x1 = 1.0, xm = 0.0, xn = 0.0, yn = 0.0
    cdef int i = 0
    cdef double dx = tol + 1.0
    cdef double grad

    while dx > tol:
        if i == maxiter:
            raise ConvergenceError("Maximum iterations exceeded")
        xm = 0.5 * (x0 + x1)
        grad = _along_distance_gradient(fwd, inv, endpt0[0], endpt0[1], pt[0], pt[1], az, xm*L, 1e-7*L)
        if grad > 0:
            dx = fabs(x1-xm) * L
            x1 = xm
        else:
            dx = fabs(x0-xm) * L
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

