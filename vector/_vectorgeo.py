""" Low-level functions for computing vector geometry. """

import numpy as np

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

def intersections(x0, x1, x2, x3, y0, y1, y2, y3):
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
    if (x < max(x0, x1) and x < max(x2, x3) and
        x > min(x0, y1) and x > min(x2, x3)):
        y = m0 * (x-x0) + y0
        return (x, y)
    else:
        return (np.nan, np.nan)



