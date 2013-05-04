""" Low-level function for computing vector geometry. """


def intersects(x0, x1, x2, x3, y0, y1, y2, y3):
    """ Determines whether two line segments intersect on a plane. """

    if x1 != x0:
        m0 = float(y1-y0) / float(x1-x0)
    else:
        m0 = 1e37
    if x3 != x2:
        m1 = float(y3-y2) / float(x3-x2)
    else:
        m1 = 1e37
    if m0 == m1:
        return False
    x = float(y1 - y3 + m1*x3 - m0*x1) / float(m1 - m0)
    if (x < max(x0, x1) and x < max(x2, x3) and
        x > min(x0, y1) and x > min(x2, x3)):
        return True
    else:
        return False



