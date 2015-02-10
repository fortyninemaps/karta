""" Low-level functions for computing vector geometry. """

import math
import numpy as np

def isbetween_inc(a, b, c):
    return b <= max(a, c) and b >= min(a, c)

def isbetween_incl(a, b, c):
    return b < max(a, c) and b >= min(a, c)

def isbetween_incr(a, b, c):
    return b <= max(a, c) and b > min(a, c)

def intersection(x0, x1, x2, x3, y0, y1, y2, y3):
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
    
    iswithinx = False
    iswithiny = False

    if abs(x1 - x0) >= 1e-15:
        if abs(x3 - x2) >= 1e-15:
            if isbetween_incl(x0, x, x1) and isbetween_incl(x2, x, x3):
                iswithinx = True
        else:
            if abs(x - x2) < 1e-15 and isbetween_incl(x0, x, x1):
                iswithinx = True
    else:
        if abs(x - x0) < 1e-15 and isbetween_incl(x2, x, x3):
            iswithinx = True

    if iswithinx:
        if abs(x-x0) >= 1e-15:
            y = m0 * (x-x0) + y0
        else:
            y = m1 * (x-x2) + y2
        if abs(y1 - y0) >= 1e-15:
            if abs(y3 - y2) >= 1e-15:
                if isbetween_incl(y0, y, y1) and isbetween_incl(y2, y, y3):
                    iswithiny = True
            elif abs(y - y2) < 1e-15 and isbetween_incl(y0, y, y1):
                iswithiny = True
        elif abs(y - y0) < 1e-15 and isbetween_incl(y2, y, y3):
            iswithiny = True
            
    if iswithinx and iswithiny:
        return (x, y)
    else:
        return (np.nan, np.nan)

def intersects_cn(xp, yp, x0, x1, y0, y1):
    """ Test whether a horizontal ray emanating left from a point (xp, yp)
    crosses a line segment. Used to implement a crossing number membership
    test.
    """
    if x1 != x0:
        m = float(y1-y0) / float(x1-x0)
    else:
        m = 1e37
    if m == 0.0:
        return False

    x = float(yp-y0) / m  + x0

    if x < xp:
        return False
    
    iswithinx = False
    iswithiny = False
    isbetween = isbetween_incl if m > 0 else isbetween_incr

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

def distance(pt0, pt1):
    """ Calculate the distance between two points (tuples) """
    d = math.sqrt(sum([abs(a-b)**2 for a, b in zip(pt0, pt1)]))
    return d

def pt_nearest_planar(pt, endpt1, endpt2):
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

def pt_nearest_proj(fwd, inv, pt, endpt0, endpt1, tol=1.0, maxiter=50):
    """ Given geodetic functions *fwd* and *inv*, a Point *pt*, and an arc from
    *endpt1* to *endpt2*, return the point on the arc that is nearest *pt*.

    Scheme employs a bisection minimization method. Iteration continues until a
    tolerance *tol* in meters is reached, or *maxiter* iterations are
    exhausted. If the iteration limit is reached, a ConvergenceError is raised.
    """
    (az, az2, L) = inv(endpt0[0], endpt0[1], endpt1[0], endpt1[1])

    def distance(x):
        trialpt = fwd(endpt0[0], endpt0[1], az, x*L)
        (_, _, d) = inv(trialpt[0], trialpt[1], pt[0], pt[1])
        return d
    
    def ddx(x):
        dx = 1e-8
        d1 = distance(x)
        d2 = distance(x+dx)
        return (d2-d1) / dx
    
    # Detect whether the nearest point is at an endpoint
    dx0, dx1 = ddx(0), ddx(1)
    if dx0 > 0:
        return endpt0, distance(0)
    elif dx1 < 0:
        return endpt1, distance(1)
    
    # Bisection iteration
    dx = tol + 1.0
    i = 0
    x0, x1 = 0, 1
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
    (xn, yn, _) = geod.fwd(endpt0[0], endpt0[1], az, xm*L)
    return (xn, yn), distance(xm)

class ConvergenceError(Exception):
    def __init__(self, message=''):
        self.message = message
    def __str__(self):
        return self.message


# QuadTree functions

def iswithin(bbox, pt):
    """ Return whether a point is within a bounding box (planar approximation). """
    if (bbox[0] <= pt[0] < bbox[1] and bbox[2] <= pt[1] < bbox[3]):
        return True
    else:
        return False

def hashpt(xmin, xmax, ymin, ymax, x, y):
    """ Returns a generator that returns successive quadrants [0-3] that
    constitute a geohash for *pt* in a global *bbox*. """
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

