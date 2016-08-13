from cpython cimport bool
from coordstring cimport CoordString
from vectorgeo cimport Vector2, azimuth, azimuth_sph, dist2
from libc.math cimport atan2

cdef bool isleft(Vector2 pt0, Vector2 pt1, Vector2 pt2):
    """ Indicate whether pt0 is to the left of a segment between pt1 and pt2.
    """
    return (pt1.x-pt0.x)*(pt2.y-pt0.y) - (pt1.y-pt0.y)*(pt2.x-pt0.x) > 0.0

def convexhull(CoordString cs):
    """ Return the convex hull for coordinates on a plane. """
    # # Find the lowermost (left?) point
    # pt0 = points[0]
    # idx = 0
    # for i, pt in enumerate(points[1:]):
    #     if (pt.y < pt0.y) or ((pt.y == pt0.y) and (pt.x < pt0.x)):
    #         pt0 = pt
    #         idx = i+1
    # points.pop(idx)

    cdef int n = len(cs)
    if cs.ring:
        n += 1
    cdef int ileftmost = -1
    cdef double x, xleftmost = cs.getX(0)
    cdef int i = 1

    while i != n:
        x = cs.getX(i)
        if x < xleftmost:
            xleftmost = x
            ileftmost = i
        elif x == xleftmost and (cs.getY(i) < cs.getY(ileftmost)):
            ileftmost = i
        i += 1

    # Sort CCW relative to pt0
    # points.sort(key=lambda pt: pt0.distance(pt))
    # points.sort(key=lambda pt: _cvectorgeo.polarangle(pt0.vertex, pt.vertex))
    # alpha = -1
    # drop = []
    # for i,pt in enumerate(points):
    #     a = _cvectorgeo.polarangle(pt0.vertex, pt.vertex)
    #     if a == alpha:
    #         drop.append(i)
    #     else:
    #         alpha = a

    # if len(drop) != 0:
    #     for i in drop[::-1]:
    #         points.pop(i)

    cdef list azimuth_indices = []
    cdef double az, y
    cdef Vector2 leftmost = Vector2(xleftmost, cs.getY(ileftmost))
    i = 0
    while i != n:
        if i == ileftmost:
            i += 1
            continue
        az = 90.0 - azimuth(leftmost, Vector2(cs.getX(i), cs.getY(i)))
        azimuth_indices.append((az, i))
        i += 1
    azimuth_indices.sort()

    # Drop all but farthest of any duplicates
    cdef int last_idx = 0
    cdef Vector2 ipt, lastpt = leftmost
    cdef list indices = [azimuth_indices[0][1]]
    i = 1
    while i != len(azimuth_indices):
        if azimuth_indices[i][0] == azimuth_indices[last_idx][0]:
            ipt = Vector2(cs.getX(azimuth_indices[i][1]), cs.getY(azimuth_indices[i][1]))
            if dist2(leftmost, ipt) > dist2(leftmost, lastpt):
                # replace last index with that of the current vertex
                indices[-1] = azimuth_indices[i][1]
                last_idx = i
                lastpt = Vector2(cs.getX(last_idx), cs.getY(last_idx))
        else:
            indices.append(azimuth_indices[i][1])
        i += 1

    # Initialize and scan for convex hull
    if len(indices) == 1:
        return [ileftmost, indices[0]]

    cdef list hull_indices = [ileftmost, indices[0], indices[1]]
    if len(indices) == 2:
        return hull_indices

    for i in indices[2:]:
        while not isleft(Vector2(cs.getX(hull_indices[-2]), cs.getY(hull_indices[-2])),
                         Vector2(cs.getX(hull_indices[-1]), cs.getY(hull_indices[-1])),
                         Vector2(cs.getX(i), cs.getY(i))):
            hull_indices.pop()
        hull_indices.append(i)
    return hull_indices

def convexhull_sph(CoordString cs):
    """ Return the convex hull for coordinates on a sphere. """
    raise NotImplementedError()

def convexhull_simple(CoordString cs):
    """ Return the convex hull for a simple polyline on a plane. """
    raise NotImplementedError()
