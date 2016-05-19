
def _flatten(vertices):
    """ Convert a nested list of coordinates into a flat list. """
    out = []
    for item in vertices:
        if hasattr(item[0], "__iter__"):
            verts = _flatten(item)
            out.extend(verts)
        else:
            out.append(item)
    return out

def _as_nested_lists(vertices):
    """ Convert a nested structure such as an ndarray into a list of lists. """
    out = []
    for part in vertices:
        if hasattr(part[0], "__iter__"):
            verts = _as_nested_lists(part)
            out.append(verts)
        else:
            out.append(list(part))
    return out

def _reproject(xy, crs1, crs2):
    """ Reproject a coordinate (or 2-tuple of x and y vectors) from *crs1* to
    *crs2*. """
    return crs1.transform(crs2, *xy)

def _reproject_nested(xy, crs1, crs2):
    """ Reproject a possibly nested list of coordinates, retaining the nesting
    structure. """
    if hasattr(xy[0], "__iter__"):
        out = []
        for xy_ in xy:
            out.append(_reproject_nested(xy_, crs1, crs2))
        return out
    else:
        return _reproject(xy, crs1, crs2)

