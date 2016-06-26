
cdef double mind(double a, double b):
    if a > b:
        return b
    else:
        return a

cdef double maxd(double a, double b):
    if a > b:
        return a
    else:
        return b

cdef int signd(double a):
    """ Return the sign of *a* """
    if a == 0.0:
        return 1
    else:
        return int(a/abs(a))

cpdef int crosses_dateline(double x0, double x1):
    """ check whether geodesic with longitudes x0, x1 [-180, 180) crosses the dateline.
    return -1 if crossing is west to east
    return 1 if crossing is east to west
    return 0 if no crossing
    """
    if (signd(x0) != signd(x1)) and (abs(x0-x1) > 180.0):
        if x1 - x0 > 180:
            return 1
        else:
            return -1
    else:
        return 0

def dateline_bbox(double[:] x, double[:] y):
    """ Return a dateline-aware bbox for geographical coordinates. """
    cdef double xmin, ymin, xmax, ymax
    cdef double rot
    cdef long i
    cdef int xdateline

    xmin = xmax = x[0]
    ymin = ymax = y[0]
    rot = 0

    for i in range(len(x)-1):
        x0, x1 = x[i], x[i+1]
        ymin = mind(ymin, y[i+1])
        ymax = maxd(ymax, y[i+1])
        xdateline = crosses_dateline(x0, x1)
        if xdateline != 0:
            rot -= xdateline * 360.0
            xmin = mind(xmin, x1+rot)
            xmax = maxd(xmax, x1+rot)
        else:
            if x0 > x1:
                xmin = mind(xmin, x1)
            else:
                xmax = maxd(xmax, x1)

    xmin = (xmin+180) % 360 - 180
    xmax = (xmax+180) % 360 - 180
    return (xmin, ymin, xmax, ymax)


