from libc.math cimport fmin, fmax, NAN
from vectorgeo cimport Vector2, fsign, bndlat_sph
from coordstring cimport CoordString

cpdef int crosses_dateline(double x0, double x1):
    """ check whether geodesic with longitudes x0, x1 [-180, 180) crosses the dateline.
    return -1 if crossing is west to east
    return 1 if crossing is east to west
    return 0 if no crossing
    """
    if (fsign(x0) != fsign(x1)) and (abs(x0-x1) > 180.0):
        if x1 - x0 > 180:
            return 1
        else:
            return -1
    else:
        return 0

def bbox(CoordString cs):
    """ Return a dateline-aware bbox for geographical coordinates. """
    cdef double xmin, ymin, xmax, ymax, segymin, segymax
    cdef double x0, y0, x1, y1
    cdef double rot
    cdef int i = 0, retint = 0, xdateline = 0, n = len(cs)

    if n == 0:
        return (NAN, NAN, NAN, NAN)

    if cs.ring:
        n += 1

    xmin = xmax = cs.getX(0)
    ymin = ymax = cs.getY(0)
    rot = 0
    x0 = cs.getX(0)
    y0 = cs.getY(0)

    while i != n:
        x1 = cs.getX(i)
        y1 = cs.getY(i)
        retint = bndlat_sph(Vector2(x0, y0), Vector2(x1, y1), &segymin, &segymax)
        if retint == 0:
            ymin = fmin(ymin, segymin)
            ymax = fmax(ymax, segymax)
        xdateline = crosses_dateline(x0, x1)
        if xdateline != 0:
            rot -= xdateline * 360.0
            xmin = fmin(xmin, x1+rot)
            xmax = fmax(xmax, x1+rot)
        else:
            if x0 > x1:
                xmin = fmin(xmin, x1)
            else:
                xmax = fmax(xmax, x1)
        x0 = x1
        y0 = y1
        i += 1

    xmin = (xmin+180) % 360 - 180
    xmax = (xmax+180) % 360 - 180
    return (xmin, ymin, xmax, ymax)
