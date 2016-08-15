from libc.math cimport fmin, fmax
from vectorgeo cimport fsign
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

def dateline_bbox(CoordString cs):
    """ Return a dateline-aware bbox for geographical coordinates. """
    cdef double xmin, ymin, xmax, ymax
    cdef double x0, x1, y1
    cdef double rot
    cdef int i
    cdef int xdateline

    xmin = xmax = cs.getX(0)
    ymin = ymax = cs.getY(0)
    rot = 0

    for i in range(len(cs)-1):
        x0 = cs.getX(i)
        x1 = cs.getX(i+1)
        y1 = cs.getY(i+1)
        ymin = fmin(ymin, y1)
        ymax = fmax(ymax, y1)
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

    xmin = (xmin+180) % 360 - 180
    xmax = (xmax+180) % 360 - 180
    return (xmin, ymin, xmax, ymax)
