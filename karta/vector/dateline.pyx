
cdef double _mind(double a, double b):
    if a > b:
        return b
    else:
        return a

cdef double _maxd(double a, double b):
    if a > b:
        return a
    else:
        return b

cdef int _signd(double a):
    """ Return the sign of *a* """
    if a == 0.0:
        return 1
    else:
        return int(a/abs(a))

cdef int _seg_crosses_dateline(double x0, double x1):
    return (_signd(x0) != _signd(x1)) and (abs(x0-x1) > 180.0)

def dateline_bbox(double[:] x, double[:] y):
    """ Return a dateline-aware bbox for geographical coordinates. """
    cdef double xmin, ymin, xmax, ymax
    cdef double rot
    cdef long i

    xmin = xmax = x[0]
    ymin = ymax = y[0]
    rot = 0

    for i in range(len(x)-1):
        x0, x1 = x[i], x[i+1]
        ymin = _mind(ymin, y[i+1])
        ymax = _maxd(ymax, y[i+1])
        if _seg_crosses_dateline(x0, x1):
            if x0 > x1:     # east to west
                rot += 360.0
            else:           # west to east
                rot -= 360.0

            xmin = _mind(xmin, x1+rot)
            xmax = _maxd(xmax, x1+rot)
        else:
            if x0 > x1:
                xmin = _mind(xmin, x1)
            else:
                xmax = _maxd(xmax, x1)

    xmin = (xmin+180) % 360 - 180
    xmax = (xmax+180) % 360 - 180
    return (xmin, ymin, xmax, ymax)


