import math
from .vector.geometry import Point
from .crs import LonLatWGS84
from collections import namedtuple

tile = namedtuple("Tile", ["z", "x", "y"])

def tile_tuple(point, zoom):
    """ Return the (z, x, y) locator for an OpenStreetMap tile containing a
    point.

    Parameters
    ----------
    point : Point
        geographic point to be contained within tile
    zoom : int
        non-negative zoom level (typically 0-18)
    """
    z = int(zoom)
    dlon = 256
    dlat = 256

    lon0, lat0 = point.crs.project(*point.vertex[:2], inverse=True)
    c = 128/math.pi * 2**z
    x0 = c * (lon0*math.pi/180+math.pi)
    y0 = c * (math.pi-math.log(math.tan(math.pi/4+lat0*math.pi/360)))

    x = int(x0 // dlon)
    y = int(y0 // dlat)
    return tile(z, x, y)

def tile_nw_corner(z, x, y):
    """ Return the northwest corner point from the tile specified by a`
    tuple.
    """
    n = 2**z
    lon = float(x)/n * 360.0 - 180.0
    lat = math.atan(math.sinh(math.pi * (1-2*y/n))) * 180.0 / math.pi
    return Point((lon, lat), crs=LonLatWGS84)

def tile_bbox(z, x, y):
    """ Return a tuple representing the bounding box of a tile tuple.
    """
    pts = [tile_nw_corner(z, x, y), tile_nw_corner(z, x+1, y),
           tile_nw_corner(z, x, y+1)]
    return pts[0].x, pts[2].y, pts[1].x, pts[0].y