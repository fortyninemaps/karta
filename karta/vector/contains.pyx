from libc.math cimport fabs
from cpython cimport bool
from coordstring cimport CoordString
from vectorgeo cimport Vector2

cdef double isleft(Vector2 pt, Vector2 pt0, Vector2 pt1) nogil:
    """ tests whether *pt* is left of segment (pt0, pt1).

    returns > 0 if left, < 0 if right, and = 0 if pt is on the segment
    """
    return (pt1.x - pt0.x) * (pt.y - pt0.y) - (pt.x - pt0.x) * (pt1.y - pt0.y)

def contains(double x, double y, CoordString poly):
    """ Uses a winding number scheme to compute whether *poly* contains a point
    (x, y) """
    if not poly.ring:
        raise TypeError("contains requires a closed CoordString")

    cdef int cnt = 0
    cdef int i
    cdef Vector2 pt, pt0, pt1

    pt = Vector2(x, y)
    pt0 = Vector2(poly.coords[0], poly.coords[1])
    for i in range(poly.length):
        if i != poly.length-1:
            pt1.x = poly.coords[(i+1)*poly.rank]
            pt1.y = poly.coords[(i+1)*poly.rank+1]
        else:
            pt1.x = poly.coords[0]
            pt1.y = poly.coords[1]
        if (pt0.y <= pt.y < pt1.y):
            # crosses upward
            if isleft(pt, pt0, pt1) > 0:
                cnt += 1
        elif (pt0.y > pt.y >= pt1.y):
            # crosses downward
            if isleft(pt, pt0, pt1) < 0:
                cnt -= 1
        pt0 = pt1

    return cnt != 0

def contains_proj(double x, double y, CoordString poly, object crs):
    """ contains implementation for geographical coordinates.
    calls crs.inverse n times, making this relatively inefficient.
    """
    cdef double sum_az = 0.0, az0, az1
    cdef double x0, y0, x1, y1
    cdef int i

    if not poly.ring:
        raise TypeError("no membership test for non-ring coordinate strings")

    x0 = poly.getX(0)
    y0 = poly.getY(0)
    az0, _, _ = crs.inverse(x, y, x0, y0)
    az0 = (az0 + 360 ) % 360
    for i in range(len(poly)):
        x1 = poly.getX(i+1)
        y1 = poly.getY(i+1)

        az1, _, _ = crs.inverse(x, y, x1, y1)
        az1 = (az1 + 360 ) % 360

        # correct for crossing prime meridian
        if az1 - az0 < -180:
            # ccw
            sum_az += (az1-az0) + 360
        elif az0 - az1 > 180:
            # cc
            sum_az += (az1-az0) - 360
        else:
            sum_az += (az1-az0)

        az0 = az1
        x0 = x1
        y0 = y1

    if fabs(sum_az) > 1e-4:
        return False
    else:
        return True

