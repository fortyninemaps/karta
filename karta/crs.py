# -*- coding: utf-8 -*-
""" Coordinate reference systems

Implements CRS classes for different kinds of spatial reference systems:

    - CartesianCRS
    - GeographicalCRS
    - SphericalCRS
    - EllipsoidalCRS
    - Proj4CRS

Coordinate reference system (CRS) objects represent mappings geographical
positions to an x, y representation, and contain both projection and geodetic
information. CRS objects should be treated is immutable. It is possible for
multiple CRS objects to implement the same system, and tests for equality may
fail.
"""

import numpy as np
import pyproj
from . import geodesy
from . import errors

try:
    import osgeo.osr
    osgeo.osr.UseExceptions()
    HASOSR = True
except ImportError:
    HASOSR = False

class CRS(object):
    """ Base class for coordinate system instances.

    Subclasses should define:

    - name attribute
    - project(x, y) method
    - forward(x, y, azimuth, distance) method
    - inverse(x0, y0, x1, y1) method

    *project* transforms between world and projected coordinates, and takes the
    optional boolean parameter *inverse* (default `False`)

    *forward* performs a forward geodetic calculation and returns world
    coordinates

    *inverse* performs an inverse geodetic calculation and returns an azimuth,
    a back azimuth, and a distance

    A CRS subclass may optionally also provide string attributes:

    - ref_proj4
    - ref_wkt

    which are used to provide interoperability with other systems.
    """
    def __str__(self):
        return "<CRS {0}>".format(self.name)

    def get_proj4(self):
        if hasattr(self, "ref_proj4"):
            return self.ref_proj4
        else:
            if not HASOSR:
                raise errors.CRSError("no ref_proj4 attribute and no conversion possible (osgeo.osr not installed)")
            srs = osgeo.osr.SpatialReference()
            if hasattr(self, "ref_wkt"):
                srs.ImportFromWkt(self.ref_wkt)
            else:
                raise AttributeError("incomplete CRS definition (missing ref_proj4, ref_wkt attributes)")
            return srs.ExportToProj4()

    def get_wkt(self):
        if hasattr(self, "ref_wkt"):
            return self.ref_proj4
        else:
            if not HASOSR:
                raise errors.CRSError("no ref_wkt attribute and no conversion possible (osgeo.osr not installed)")
            srs = osgeo.osr.SpatialReference()
            if hasattr(self, "ref_proj4"):
                srs.ImportFromProj4(self.ref_proj4.replace("lon", "long"))
            else:
                raise AttributeError("incomplete CRS definition (missing ref_proj4, ref_wkt attributes)")
            return srs.ExportToWkt()

class CartesianCRS(CRS):
    """ Cartesian (flat-earth) reference systems with (x, y) coordinates """
    name = "Cartesian"

    ref_proj4 = ""
    ref_wkt = ""

    @staticmethod
    def project(x, y, inverse=False):
        # Projection on a cartesian coordinate system is the identity
        return x, y

    @staticmethod
    def forward(x, y, az, dist, radians=False):
        """ Forward geodetic problem from a point """
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
        """ Inverse geodetic problem to find the geodesic between points """
        dist = geodesy.plane_distance(x1, y1, x2, y2)
        az = geodesy.plane_azimuth(x1, y1, x2, y2)
        baz = geodesy.unroll_angle(az + np.pi)

        if not radians:
            az = az * 180 / np.pi
            baz = baz * 180 / np.pi
        return az, baz, dist

class GeographicalCRS(CRS):
    """ Reference systems with longitude-latitude (θ, φ) coordinates.

    `spheroid` refers to a proj.4 spheroid identifier, e.g. "+ellps=WGS84"
    """
    def __init__(self, spheroid, name):
        self._geod = pyproj.Geod(spheroid)
        self.name = name
        self.ref_proj4 = "+proj=lonlat %s" % self._geod.initstring
        return

    @staticmethod
    def project(x, y, inverse=False, radians=False):
        # Projection on a geographical coordinate system is the identity
        if not radians:
            return x, y
        else:
            return np.array(x)/180 * np.pi, np.array(y)/180 * np.pi

    def forward(self, *args, **kwargs):
        """ Forward geodetic problem from a point """
        return self._geod.fwd(*args, **kwargs)

    def inverse(self, *args, **kwargs):
        """ Inverse geodetic problem to find the geodesic between points """
        return self._geod.inv(*args, **kwargs)

class SphericalCRS(GeographicalCRS):
    """ Spherical geographic coordinate system defined by a radius. """
    name = "Spherical"

    def __init__(self, radius):
        self.radius = radius
        return

    @property
    def ref_proj4(self):
        return "+proj=lonlat +a=%f +b=%f +no_defs" % (self.radius, self.radius)

    @property
    def ref_wkt(self):
        return ('GEOGCS["unnamed ellipse",DATUM["unknown",'
                'SPHEROID["unnamed",%f,0]],PRIMEM["Greenwich",0],'
                'UNIT["degree",0.0174532925199433]]' % self.radius)

    def forward(self, lons, lats, az, dist, radians=False):
        """ Forward geodetic problem from a point """
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
        """ Inverse geodetic problem to find the geodesic between points """
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

class EllipsoidalCRS(GeographicalCRS):
    """ Ellipsoidal geographic coordinate system defined by equatorial and
    polar radii.
    """

    def __init__(self, a, b):
        """ Define a geographical coordinate system on an ellipse. *a* is the
        equatorial radius, and *b* is the polar radius. """
        self.a = a
        self.b = b
        self.ref_proj4 = "+proj=lonlat +units=m +a=%f +b=%f +no_defs" % (a, b)
        return

    def get_proj4(self):
        return self.ref_proj4

    def project(self, x, y, **kw):
        # Projection on a geographical coordinate system is the identity
        return x, y

    def forward(self, x, y, azimuth, distance, radians=False):
        """ Forward geodetic problem from a point """
        if radians:
            x *= 180.0/pi
            y *= 180.0/pi
            azimuth *= 180.0/pi

        x2, y2, baz = geodesy.ellipsoidal_forward(self.a, self.b, x, y, azimuth, distance)

        if radians:
            x2 *= pi/180.0
            y2 *= pi/180.0
            baz *= pi/180.0

        return x2, y2, baz

    def inverse(self, x1, y1, x2, y2, radians=False):
        """ Inverse geodetic problem to find the geodesic between points """
        if radians:
            x1 *= 180.0/pi
            y1 *= 180.0/pi
            x2 *= 180.0/pi
            y2 *= 180.0/pi

        az, baz, dist = geodesy.ellipsoidal_inverse(self.a, self.b, x1, y1, x2, y2)

        if radians:
            az *= pi/180.0
            baz *= pi/180.0

        return az, baz, dist

class Proj4CRS(CRS):
    """ Custom reference systems, which may be backed by a *pypoj.Proj* instance
    or a custom projection function.

    `proj` is a dictionary used to define an instance of `pyproj.Proj`

    `spheroid` refers to a proj.4 spheroid identifier, e.g. "+ellps=WGS84"
    """
    def __init__(self, proj, spheroid=None, name=None):
        if spheroid is None:
            # Extract the spheroid from the +ellps, +datum, or +a,+b parameters
            # if present
            if "+ellps" in proj:
                for kv in proj.split():
                    if "+ellps" in kv:
                        k,v = kv.split("=")
                        spheroid = "+ellps=%s" % v
                        break
            elif "+datum" in proj:
                for kv in proj.split():
                    if "+datum" in kv:
                        k,v = kv.split("=")
                        spheroid = "+ellps=%s" % v
                        break
            elif ("+a" in proj) and ("+b" in proj):
                a, b = None, None
                for kv in proj.split():
                    if "+a" in kv:
                        _,a = kv.split("=")
                    elif "+b" in kv:
                        _,b = kv.split("=")
                    if None not in (a,b):
                        spheroid = "+a=%s +b=%s" % (a, b)
                        break
            else:
                raise errors.CRSError("Spheroid must be provided: %s" % proj)

        self.project = pyproj.Proj(proj)
        self._geod = pyproj.Geod(spheroid)

        self.initstring_proj = self.project.srs
        self.initstring_geod = self._geod.initstring
        self.ref_proj4 = "%s %s" % (self.project.srs, self._geod.initstring)

        if name is not None:
            self.name = name
        else:
            self.name = proj
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
        """ Forward geodetic problem from a point """
        return self._geod.fwd(*args, **kwargs)

    def inverse(self, *args, **kwargs):
        """ Inverse geodetic problem to find the geodesic between points """
        return self._geod.inv(*args, **kwargs)

def crs_from_wkt(wkt):
    srs = osgeo.osr.SpatialReference(wkt)
    srs.MorphFromESRI()     # munge parameters for proj.4 export
    if srs.IsGeographic():
        proj4 = srs.ExportToProj4()
        for arg in proj4.split():
            if "+ellps" in arg:
                return GeographicalCRS(arg, "unnamed")
            if "+datum" in arg:
                _,v = arg.split("=")
                return GeographicalCRS("+ellps=%s" % v, "unnamed")
        raise errors.CRSError("Could not interpret %s as geographical CRS" % proj4)
    else:
        return Proj4CRS(srs.ExportToProj4())

############ Predefined CRS instances ############

Cartesian = CartesianCRS()
SphericalEarth = SphericalCRS(6371009.0)
LonLatWGS84 = EllipsoidalCRS(6378137.0, 6356752.314245)
LonLatNAD83 = EllipsoidalCRS(6378137.0, 6356752.314140)
LonLatNAD27 = EllipsoidalCRS(6378206.4, 6356583.8)

LonLatWGS84_proj4 = GeographicalCRS("+ellps=WGS84", "WGS84 (Geographical)")
LonLatNAD27_proj4 = GeographicalCRS("+ellps=clrk66", "NAD27 (Geographical)")
LonLatNAD83_proj4 = GeographicalCRS("+ellps=GRS80", "NAD83 (Geographical)")

UPSNorth = Proj4CRS(proj="+proj=stere +lat_0=90 +lat_ts=90 +lon_0=0 +k=0.994 +x_0=2000000 +y_0=2000000 +units=m +datum=WGS84 +no_defs",
        spheroid="+ellps=WGS84", name="Universal Polar Stereographic (North)")

UPSSouth = Proj4CRS(proj="+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 +k=0.994 +x_0=2000000 +y_0=2000000 +units=m +datum=WGS84 +no_defs",
        spheroid="+ellps=WGS84", name="Universal Polar Stereographic (South)")

NSIDCNorth = Proj4CRS(proj="+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +units=m +datum=WGS84 +no_defs",
        spheroid="+ellps=WGS84", name="NSIDC (North)")

NSIDCSouth = Proj4CRS(proj="+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 +k=1 +x_0=0 +y_0=0 +units=m +datum=WGS84 +no_defs",
        spheroid="+ellps=WGS84", name="NSIDC (South)")

LambertEqualArea = Proj4CRS(proj="+proj=laea +lat_0=0 +lon_0=0 +x_0=0 +y_0=0",
        spheroid="+ellps=WGS84", name="Lambert Equal Area")

GallPetersEqualArea = Proj4CRS("+proj=cea +lon_0=0 +lat_ts=45 +x_0=0 +y_0=0 "
                               "+ellps=WGS84 +units=m +no_defs",
        spheroid="+ellps=WGS84", name="Gall Peters Equal Area")

WebMercator = Proj4CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 "
                       "+lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m "
                       "+nadgrids=@null +wktext +no_defs",
        spheroid="+a=6378137 +b=6378137", name="Web Mercator")

