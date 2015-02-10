# -*- coding: utf-8 -*-
"""
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
from . import geodesy

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

class Cartesian(CRS):
    name = "Cartesian (flat-earth)"

    @staticmethod
    def project(x, y, inverse=False):
        return x, y

    @staticmethod
    def forward(x, y, az, dist, radians=False):
        """ Returns x, y, and back azimuths """
        if not radians:
            az = np.array(az) / 180 * np.pi

        x2 = x + dist * np.sin(az)
        y2 = y + dist * np.cos(az)
        baz = geodesy.unroll_angle(az + np.pi)

        if not radians:
            baz = np.array(baz) * 180 / np.pi
        return x2, y2, baz

    @staticmethod
    def inverse(x1, y1, x2, y2, radians=False):
        """ Returns forward and back azimuths and distances """
        dist = geodesy.plane_distance(x1, y1, x2, y2)
        az = geodesy.plane_azimuth(x1, y1, x2, y2)
        baz = geodesy.unroll_angle(az + np.pi)

        if not radians:
            az = az * 180 / np.pi
            baz = baz * 180 / np.pi
        return az, baz, dist

class Spherical(CRS):
    """ Spherical (θ, φ) coordinate system. """
    name = "Spherical"

    def __init__(self, radius):
        self.radius = radius

    @staticmethod
    def project(x, y, inverse=False, radians=False):
        if radians:
            return x, y
        else:
            return np.array(x)/180 * np.pi, np.array(y)/180 * np.pi

    def forward(self, lons, lats, az, dist, radians=False):
        """ Returns lons, lats, and back azimuths """
        if not radians:
            lons = np.array(lons) * np.pi / 180.0
            lats = np.array(lats) * np.pi / 180.0
            az = np.array(az) * np.pi / 180.0


        d_ = dist / self.radius
        lats2 = np.arcsin(np.sin(lats) * np.cos(d_) + 
                    np.cos(lats) * np.sin(d_) * np.cos(az))
        dlons = np.arccos((np.cos(d_) - np.sin(lats2) * np.sin(lats)) /
                          (np.cos(lats) * np.cos(lats2)))
        baz = np.arccos((np.sin(lats) - np.cos(d_) * np.sin(lats2)) /
                        (np.sin(d_) * np.cos(lats2)))
        if 0 <= az < np.pi:
            lons2 = lons + dlons
            baz = -baz
        elif np.pi <= az < 2*np.pi:
            lons2 = lons - dlons
        else:
            raise ValueError("azimuth should be [0, 2pi)")

        baz = geodesy.unroll_angle(baz)

        if not radians:
            lons2 = np.array(lons2) * 180 / np.pi
            lats2 = np.array(lats2) * 180 / np.pi
            baz = np.array(baz) * 180 / np.pi
        return lons2, lats2, baz

    def inverse(self, lons1, lats1, lons2, lats2, radians=False):
        """ Returns forward and back azimuths and distances """
        if not radians:
            lons1 = lons1 * np.pi / 180.0
            lons2 = lons2 * np.pi / 180.0
            lats1 = lats1 * np.pi / 180.0
            lats2 = lats2 * np.pi / 180.0

        az = geodesy.sphere_azimuth(lons1, lats1, lons2, lats2)
        baz = geodesy.sphere_azimuth(lons2, lats2, lons1, lats1)
        dist = geodesy.sphere_distance(lons1, lats1, lons2, lats2, self.radius)

        if not radians:
            az = az * 180 / np.pi
            baz = baz * 180 / np.pi
        return az, baz, dist

SphericalEarth = Spherical(6371.0)

class CustomCRS(CRS):

    def __init__(self, epsg=None, proj=None, geod=None, wkt=None, urn=None, name=None):
        if epsg is not None:
            raise NotImplementedError("EPSG lookup not implemented")
        elif None not in (proj, geod):
            self.project = pyproj.Proj(proj)
            self._geod = pyproj.Geod(geod)
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
        if (getattr(self.project, "srs", 0) == getattr(other.project, "srs", 1) and
            getattr(self._geod, "initstring", 0) == getattr(other._geod, "initstring", 1)):
            return True
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def forward(self, *args, **kwargs):
        return self._geod.fwd(*args, **kwargs)

    def inverse(self, *args, **kwargs):
        return self._geod.inv(*args, **kwargs)

class CRSError(Exception):
    """ Exception to raise for invalid geodetic operations. """
    def __init__(self, message=''):
        self.message = message

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


