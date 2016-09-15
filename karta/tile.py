import math
from .vector.geometry import Point
from .crs import LonLatWGS84

class Tile(object):

    def __init__(self, z, x, y):
        self.z = z
        self.x = x
        self.y = y

    def __eq__(self, other):
        if not isinstance(other, Tile):
            return False
        return (self.z == other.z) and (self.x == other.x) and (self.y == other.y)
    
    def __neq__(self, other):
        return not (self == other)

    def nw_corner(self):
        """ Return the northwest corner point from a Tile.

        Returns
        -------
        Point
        """
        z = self.z
        x = self.x
        y = self.y
        n = 2**z
        lon = float(x)/n * 360.0 - 180.0
        lat = math.atan(math.sinh(math.pi * (1-2*y/n))) * 180.0 / math.pi
        return Point((lon, lat), crs=LonLatWGS84)

    @property
    def bbox(self):
        """ Return a tuple representing the bounding box of a Tile

        Returns
        -------
        tuple
            (xmin, ymin, xmax, ymax)
        """
        z = self.z
        x = self.x
        y = self.y
        pts = [self.nw_corner(),
               Tile(z, x+1, y).nw_corner(),
               Tile(z, x, y+1).nw_corner()]
        return pts[0].x, pts[2].y, pts[1].x, pts[0].y

def tile_from_point(point, zoom):
    """ Return the (z, x, y) locator for an OpenStreetMap tile containing a
    point.

    Parameters
    ----------
    point : Point
        geographic point to be contained within tile
    zoom : int
        non-negative zoom level (typically 0-18)

    Returns
    -------
    Tile
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
    return Tile(z, x, y)
