# -*- coding: utf-8 -*-
"""
Experimental re-write of karta's CRS system.

Desired interface:

    from karta.crs import Cartesian
    from karta.crs import Spherical
    from karta.crs import LonLatWGS84

    from karta.crs import CustomCRS

    crs = CustomCRS(wkt=\"\"\"PROJCS["WGS 84 / NSIDC Sea Ice Polar Stereographic North",
    GEOGCS["WGS 84",
        DATUM["WGS_1984",
            SPHEROID["WGS 84",6378137,298.257223563,
                AUTHORITY["EPSG","7030"]],
            AUTHORITY["EPSG","6326"]],
        PRIMEM["Greenwich",0,
            AUTHORITY["EPSG","8901"]],
        UNIT["degree",0.01745329251994328,
            AUTHORITY["EPSG","9122"]],
        AUTHORITY["EPSG","4326"]],
    UNIT["metre",1,
        AUTHORITY["EPSG","9001"]],
    PROJECTION["Polar_Stereographic"],
    PARAMETER["latitude_of_origin",70],
    PARAMETER["central_meridian",-45],
    PARAMETER["scale_factor",1],
    PARAMETER["false_easting",0],
    PARAMETER["false_northing",0],
    AUTHORITY["EPSG","3413"],
    AXIS["X",UNKNOWN],
    AXIS["Y",UNKNOWN]]\"\"\")
    
    crs = CustomCRS(proj4=\"\"\"+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs\"\"\")

    crs = CustomCRS(epsg=3413)
"""

import numpy as np
import pyproj

# A CRS class needs to have pyproj Proj and Geod instances. Exceptions are
# `Cartesian` and `Spherical` which are implemented specially.

class CRS(object):
    """ Base class for coordinate system instances
    
    Subclasses should at a minimum define a *name* attribute. Providing *proj*
    and *geod* methods permits a geometry associated with a CRS subclass to be
    used in analyses.
    """
    def __str__(self):
        return "<CRS {0}>".format(self.name)

class _cartesian_geod(object):
    """ Static class that substitutes for a pyrproj.Geod instance for
    cartesian coordinates """

    @staticmethod
    def fwd(lons, lats, az, dist, radians=False):
        """ Returns lons, lats, and back azimuths """
        if radians:
            lons = np.array(lons) * 180 / np.pi
            lats = np.array(lats) * 180 / np.pi
        else:
            az = np.array(az) / 180 * np.pi

        backaz = az + np.pi
        if hasattr(backaz, "__len__"):
            m = backaz >= 2*np.pi
            backaz[m] -= 2*np.pi
        elif backaz >= 2*np.pi:
            backaz -= 2*np.pi

        lons2 = lons + dist * np.sin(az)
        lats2 = lats + dist * np.cos(az)

        if radians:
            lons2 = np.array(lons2) / 180 * np.pi
            lats2 = np.array(lats2) / 180 * np.pi
        else:
            backaz = np.array(backaz) * 180 / np.pi
        return lons2, lats2, backaz

    @staticmethod
    def inv(lons1, lats1, lons2, lats2, radians=False):
        """ Returns forward and back azimuths and distances """
        dy = np.array(lats2) - np.array(lats1)
        dx = np.array(lons2) - np.array(lons1)

        dist = np.sqrt(dx**2 + dy**2)
        az = plane_azimuth(dx, dy)
        backaz = az + np.pi
        if hasattr(backaz, "__len__"):
            m = backaz >= 2*np.pi
            backaz[m] -= 2*np.pi
        elif backaz >= 2*np.pi:
            backaz -= 2*np.pi

        if not radians:
            az = az * 180 / np.pi
            backaz = backaz * 180 / np.pi
        return az, backaz, dist

class _spherical_geod(object):
    """ Static class that substitutes for a pyrproj.Geod instance for
    spherical coordinates """

    def __init__(self, radius):
        self.radius = radius

    @staticmethod
    def _distance(lons1, lats1, lons2, lats2):
        """ Computes the great circle distance between point pairs on a
        sphere. 

        REQUIRES RADIANS
        """
        if dx > 0.01 or dy > 0.01:
            # use spherical law of cosines
            d_ = np.arccos(np.sin(lats1) * np.sin(lats2) +
                    np.cos(lats1) * np.cos(lats2) * np.cos(np.abs(lons2-lons1)))

        else:
            # use haversine
            dy = np.abs(lats2 - lats1)
            dx = np.abs(lons2 - lons1)

            d_ = (2 * np.arcsin(
                        np.sqrt(np.sin(dy / 2.)**2 +
                            np.cos(lats1) *
                            np.cos(lats2) *
                            np.sin(dx / 2.)**2)))
        return self.radius * d_

    @staticmethod
    def fwd(lons, lats, az, dist, radians=False):
        """ Returns lons, lats, and back azimuths """
        if not radians:
            lons = np.array(lons) * np.pi / 180.0
            lats = np.array(lats) * np.pi / 180.0
            az = np.array(az) * np.pi / 180.0


        d_ = dist / self.radius
        lats2 = np.arccos(np.sin(lats) * np.cos(d_) + 
                    np.cos(lats) * np.sin(d_) * np.sin(az))
        lons2 = lons + np.arccos((np.cos(d_) - np.sin(lats2) * np.sin(lats)) /
                                 np.cos(lats) * np.cos(lats2))
        backaz = az + np.pi
        if hasattr(backaz, "__len__"):
            m = backaz >= 2*np.pi
            backaz[m] -= 2*np.pi
        elif backaz >= 2*np.pi:
            backaz -= 2*np.pi

        if not radians:
            lons2 = np.array(lons2) * 180 / np.pi
            lats2 = np.array(lats2) * 180 / np.pi
            backaz = np.array(backaz) * 180 / np.pi
        return lons2, lats2, backaz

    @staticmethod
    def inv(lons1, lats1, lons2, lats2, radians=False):
        """ Returns forward and back azimuths and distances """
        if not radians:
            lons1 = lons1 * np.pi / 180.0
            lons2 = lons2 * np.pi / 180.0
            lats1 = lats1 * np.pi / 180.0
            lats2 = lats2 * np.pi / 180.0

        dlon = lons2-lons1
        az = np.arctan(np.sin(dlon) / 
                (np.cos(lats1) * np.tan(lats2) - np.sin(lats1) * np.cos(dlon)))

        backaz = az + np.pi
        if hasattr(backaz, "__len__"):
            m = backaz >= 2*np.pi
            backaz[m] -= 2*np.pi
        elif backaz >= 2*np.pi:
            backaz -= 2*np.pi
            
        dist = self._distance(lons1, lats1, lons2, lats2)

        if not radians:
            az = az * 180 / np.pi
            backaz = backaz * 180 / np.pi
        return az, backaz, dist

class Cartesian(CRS):
    name = "Cartesian (flat-earth)"
    geod = _cartesian_geod

    @staticmethod
    def proj(x, y, inverse=False):
        return x, y

class Spherical(CRS):
    """ Spherical (θ, φ) coordinate system. """
    _geod = _spherical_geod
    def __init__(self, radius):
        self.radius = radius

    @staticmethod
    def proj(x, y, inverse=False, radians=False):
        if radians:
            return x, y
        else:
            return np.array(x)/180 * np.pi, np.array(y)/180 * np.pi

SphericalEarth = Spherical(6370.0)

class CustomCRS(CRS):

    def __init__(self, epsg=None, proj=None, geod=None, wkt=None, urn=None, name=None):
        if epsg is not None:
            raise NotImplementedError("EPSG lookup not implemented")
        elif None not in (proj, geod):
            self.proj = pyproj.Proj(proj)
            self.geod = pyproj.Geod(geod)
        elif wkt is not None:
            raise NotImplementedError("WKT lookup not implemented")
        elif urn is not None:
            raise NotImplementedError("URN lookup not implemented")

        if name is not None:
            self.name = name
        else:
            self.name = "Unnamed CRS"
        return

    def __eq__(self, other):
        if (getattr(self.proj, "srs", 0) == getattr(other.proj, "srs", 1) and
            getattr(self.geod, "initstring", 0) == getattr(other.geod, "initstring", 1)):
            return True
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

class CRSError(Exception):
    """ Exception to raise for invalid geodetic operations. """
    def __init__(self, message=''):
        self.message = message

# The following functions are here on a temporary basis
def _azimuth(dx, dy):
    """ Return cartesian azimuth between points (scalar func) """
    if dx > 0:
        if dy > 0:
            return np.arctan(dx/dy)
        elif dy < 0:
            return np.pi - np.arctan(-dx/dy)
        else:
            return 0.5*np.pi
    elif dx < 0:
        if dy > 0:
            return 2*np.pi - np.arctan(-dx/dy)
        elif dy < 0:
            return np.pi + np.arctan(dx/dy)
        else:
            return 1.5*np.pi
    else:
        if dy > 0:
            return 0.0
        else:
            return np.pi

def plane_azimuth(dx, dy):
    """ Return cartesian azimuth between points """
    if hasattr(dx, "__iter__"):
        return numpy.array([_azimuth(dx_, dy_) for dx_, dy_ in zip(dx, dy)])
    else:
        return _azimuth(dx, dy)

############ Predefined CRS instances ############

LonLatWGS84 = CustomCRS(proj="+proj=longlat +datum=WGS84 +no_defs", geod="+ellps=WGS84", name="WGS84 (Geographical)")
s = "+proj=longlat +datum=WGS84 +no_defs"
LonLatNAD27 = CustomCRS(proj="+proj=longlat +datum=WGS84 +no_defs", geod="+ellps=clrk66", name="NAD27 (Geographical)")
LonLatNAD83 = CustomCRS(proj="+proj=longlat +datum=WGS84 +no_defs", geod="+ellps=GRS80", name="NAD83 (Geographical)")

UPSNorth = CustomCRS(proj="+proj=stere +lat_0=90 +lat_ts=90 +lon_0=0 +k=0.994 +x_0=2000000 +y_0=2000000 +units=m +datum=WGS84 +no_defs",
        geod="+ellps=WGS84", name="Universal Polar Stereographic (North)")

UPSSouth = CustomCRS(proj="+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=0.994 +x_0=2000000 +y_0=2000000 +units=m +datum=WGS84 +no_defs",
        geod="+ellps=WGS84", name="Universal Polar Stereographic (South)")

NSIDCNorth = CustomCRS(proj="+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +units=m +datum=WGS84 +no_defs",
        geod="+ellps=WGS84", name="NSIDC (North)")

NSIDCSouth = CustomCRS(proj="+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +units=m +datum=WGS84 +no_defs",
        geod="+ellps=WGS84", name="NSIDC (South)")


