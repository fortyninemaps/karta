
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
    else:
        return _reproject(xy, crs1, crs2)

