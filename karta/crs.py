# -*- coding: utf-8 -*-
""" Coordinate reference systems

Implements CRS classes for different kinds of spatial reference systems:

    - CartesianCRS
    - GeographicalCRS
    - SphericalCRS
    - EllipsoidalCRS
    - ProjectedCRS

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

# proj.4 keys for ellipsoids used by named datum definitions
DATUM_ELLIPSOIDS = {"WGS84": "WGS84", "GGRS87": "GRS80", "NAD83": "GRS80",
                    "NAD27": "clrk66", "potsdam": "bessel", "carthage": "clark80",
                    "hermannskogel": "bessel", "ire65": "mod_airy",
                    "nzgd49": "intl", "OSGB36": "airy"}

ELLIPSOID_DATA = {
 "MERIT": (6378137.0, None, 298.257, "MERIT 1983"),
 "SGS85": (6378136.0, None, 298.257, "Soviet Geodetic System 85"),
 "GRS80": (6378137.0, None, 298.257222101, "GRS 1980(IUGG, 1980)"),
 "IAU76": (6378140.0, None, 298.257, "IAU 1976"),
 "airy": (6377563.396, 6356256.910, None, "Airy 1830"),
 "APL4.9": (6378137.0, None, 298.25, "Appl. Physics. 1965"),
 "NWL9D": (6378145.0, None, 298.25, "Naval Weapons Lab., 1965"),
 "mod_airy": (6377340.189, 6356034.446, None, "Modified Airy"),
 "andrae": (6377104.43, None, 300.0, "Andrae 1876 (Den., Iclnd.)"),
 "aust_SA": (6378160.0, None, 298.25, "Australian Natl & S. Amer. 1969"),
 "GRS67": (6378160.0, None, 298.2471674270, "GRS 67(IUGG 1967)"),
 "bessel": (6377397.155, None, 299.1528128, "Bessel 1841"),
 "bess_nam": (6377483.865, None, 299.1528128, "Bessel 1841 (Namibia)"),
 "clrk66": (6378206.4, 6356583.8, None, "Clarke 1866"),
 "clrk80": (6378249.145, None, 293.4663, "Clarke 1880 mod."),
 "clrk80ign": (6378249.2, None, 293.4660212936269, "Clarke 1880 (IGN)."),
 "CPM": (6375738.7, None, 334.29, "Comm. des Poids et Mesures 1799"),
 "delmbr": (6376428., None, 311.5, "Delambre 1810 (Belgium)"),
 "engelis": (6378136.05, None, 298.2566, "Engelis 1985"),
 "evrst30": (6377276.345, None, 300.8017, "Everest 1830"),
 "evrst48": (6377304.063, None, 300.8017, "Everest 1948"),
 "evrst56": (6377301.243, None, 300.8017, "Everest 1956"),
 "evrst69": (6377295.664, None, 300.8017, "Everest 1969"),
 "evrstSS": (6377298.556, None, 300.8017, "Everest (Sabah & Sarawak)"),
 "fschr60": (6378166., None, 298.3, "Fischer (Mercury Datum) 1960"),
 "fschr60m": (6378155., None, 298.3, "Modified Fischer 1960"),
 "fschr68": (6378150., None, 298.3, "Fischer 1968"),
 "helmert": (6378200., None, 298.3, "Helmert 1906"),
 "hough": (6378270.0, None, 297., "Hough"),
 "intl": (6378388.0, None, 297., "International 1909 (Hayford)"),
 "krass": (6378245.0, None, 298.3, "Krassovsky, 1942"),
 "kaula": (6378163., None, 298.24, "Kaula 1961"),
 "lerch": (6378139., None, 298.257, "Lerch 1979"),
 "mprts": (6397300., None, 191., "Maupertius 1738"),
 "new_intl": (6378157.5, 6356772.2, None, "New International 1967"),
 "plessis": (6376523., 6355863., None, "Plessis 1817 (France)"),
 "SEasia": (6378155.0, 6356773.3205, None, "Southeast Asia"),
 "walbeck": (6376896.0, 6355834.8467, None, "Walbeck"),
 "WGS60": (6378165.0, None, 298.3, "WGS 60"),
 "WGS66": (6378145.0, None, 298.25, "WGS 66"),
 "WGS72": (6378135.0, None, 298.26, "WGS 72"),
 "WGS84": (6378137.0, None, 298.257223563, "WGS 84"),
 "sphere": (6370997.0, 6370997.0, None, "Normal Sphere (r=6370997)")}

def ellipsoid_major_axis(name):
    a, _, _ = ELLIPSOID_DATA[name]
    return a

def ellipsoid_minor_axis(name):
    a, b, rf = ELLIPSOID_DATA[name]
    if b is None:
        b = a-a/rf
    return b

def ellipsoid_flattening(name):
    a, b, rf = ELLIPSOID_DATA[name]
    if rf is None:
        rf = a/(a-b)
    return 1.0/rf

class Ellipsoid(object):
    def __init__(self, name, a=None, b=None, f=None, rf=None):
        if a is None:
            raise ValueError("major axis (a) must be provided")

        if b is not None:
            f = (a-b)/a
        elif f is not None:
            b = a-a*f
        elif rf is not None:
            f = 1.0/rf
            b = a-a*f

        self.name = name
        self.a = a
        self.b = b
        self.f = f
        return

class CRS(object):
    """ Base class for coordinate system instances.

    Subclasses should define:

    - name attribute
    - project(x, y) method
    - transform(other, x, y) method
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
                raise RuntimeError("no ref_proj4 attribute and no conversion possible (osgeo.osr not installed)")
            srs = osgeo.osr.SpatialReference()
            if hasattr(self, "ref_wkt"):
                srs.ImportFromWkt(self.ref_wkt)
            else:
                raise AttributeError("incomplete CRS definition (missing ref_proj4, ref_wkt attributes)")
            return srs.ExportToProj4()

    def get_wkt(self):
        if hasattr(self, "ref_wkt"):
            return self.ref_wkt
        else:
            if not HASOSR:
                raise RuntimeError("no ref_wkt attribute and no conversion possible (osgeo.osr not installed)")
            srs = osgeo.osr.SpatialReference()
            if hasattr(self, "ref_proj4"):
                _proj4 = self.ref_proj4.replace("latlon", "latlong")\
                                       .replace("lonlat", "longlat")
                srs.ImportFromProj4(_proj4)
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
        ela, elb, elrf, ename = ELLIPSOID_DATA[spheroid.split("=")[1]]
        self.ellipsoid = Ellipsoid(ename, ela, b=elb, rf=elrf)
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
        x, y, baz = self._geod.fwd(*args, **kwargs)
        baz = (baz + 180) % 360 - 180
        return x, y, baz

    def inverse(self, *args, **kwargs):
        """ Inverse geodetic problem to find the geodesic between points """
        az, baz, dist = self._geod.inv(*args, **kwargs)
        az = (az + 180) % 360 - 180
        baz = (baz + 180) % 360 - 180
        return az, baz, dist

class ProjectedCRS(CRS):
    """ Custom reference systems, which may be backed by a *pypoj.Proj* instance
    or a custom projection function.

    `proj` is a dictionary or string used to define an instance of `pyproj.Proj`

    `spheroid` refers to a proj.4 spheroid identifier, e.g. "+ellps=WGS84"
    """
    def __init__(self, proj, spheroid=None, name=None):
        if spheroid is None:
            ellipsoid = parse_ellipsoid(proj)
        else:
            ellipsoid = parse_ellipsoid(spheroid)

        self.project = pyproj.Proj(proj)
        self._geod = pyproj.Geod("+a=%s +b=%s" % (ellipsoid.a, ellipsoid.b))
        self.ellipsoid = ellipsoid

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
        x, y, baz = self._geod.fwd(*args, **kwargs)
        baz = (baz + 180) % 360 - 180
        return x, y, baz

    def inverse(self, *args, **kwargs):
        """ Inverse geodetic problem to find the geodesic between points """
        az, baz, dist = self._geod.inv(*args, **kwargs)
        az = (az + 180) % 360 - 180
        baz = (baz + 180) % 360 - 180
        return az, baz, dist

    def transform(self, other, x, y, z=None):
        return pyproj.transform(self.project, other.project, x, y, z=z)

class SphericalCRS(GeographicalCRS):
    """ Spherical geographic coordinate system defined by a radius.
    
    DEPRECATED 
    """
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

    DEPRECATED
    """
    name = "Ellipsoidal"

    def __init__(self, a, b):
        """ Define a geographical coordinate system on an ellipse. *a* is the
        equatorial radius, and *b* is the polar radius. """
        self.a = a
        self.b = b
        self.ref_proj4 = "+proj=lonlat +units=m +a=%f +b=%f +no_defs" % (a, b)
        return

    def forward(self, x, y, azimuth, distance, radians=False):
        """ Forward geodetic problem from a point """
        if radians:
            x *= 180.0/pi
            y *= 180.0/pi
            azimuth *= 180.0/pi

        x2, y2, baz = geodesy.ellipsoidal_forward(self.a, self.b, x, y, azimuth, distance)
        baz = (baz + 180) % 360 - 180

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
        az = (az + 180) % 360 - 180
        baz = (baz + 180) % 360 - 180

        if radians:
            az *= pi/180.0
            baz *= pi/180.0

        return az, baz, dist

def parse_ellipsoid(projstring):
    ename, ela, elb, elrf = None, None, None, None
    if "+ellps" in projstring:
        for kv in projstring.split():
            if kv.startswith("+ellps"):
                k,v = kv.split("=")
                ela, elb, elrf, ename = ELLIPSOID_DATA[v]
                break
    elif "+datum" in projstring:
        for kv in projstring.split():
            if kv.startswith("+datum"):
                k,v = kv.split("=")
                ellps = DATUM_ELLIPSOIDS[v]
                ela, elb, elrf, ename = ELLIPSOID_DATA[ellps]
                break
    elif ("+a" in projstring) and (("+b" in projstring) or ("+f" in projstring)):
        ename = "Unknown"
        for kv in projstring.split():
            if kv.startswith("+a"):
                _, ela = kv.split("=")
                ela = float(ela)
            elif kv.startswith("+b"):
                _, elb = kv.split("=")
                elb = float(elb)
            elif kv.startswith("+f"):
                _, elf = kv.split("=")
                elrf = 1.0/float(elf)
            if ela and (elb or elrf):
                break
    else:
        raise errors.CRSError("ellipsoid could not be extracted from %s" % projstring)
    return Ellipsoid(ename, ela, b=elb, rf=elrf)


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
                return GeographicalCRS("+ellps=%s" % DATUM_ELLIPSOIDS[v], v)
        raise errors.CRSError("Could not interpret %s as geographical CRS" % proj4)
    else:
        return ProjectedCRS(srs.ExportToProj4())

############ Predefined CRS instances ############

Cartesian = CartesianCRS()

_SphericalEarth = SphericalCRS(6371009.0)
_LonLatWGS84 = EllipsoidalCRS(6378137.0, 6356752.314245)
_LonLatNAD83 = EllipsoidalCRS(6378137.0, 6356752.314140)
_LonLatNAD27 = EllipsoidalCRS(6378206.4, 6356583.8)

SphericalEarth = GeographicalCRS("+ellps=sphere", "Normal Sphere")
LonLatWGS84 = GeographicalCRS("+ellps=WGS84", "WGS84 (Geographical)")
LonLatNAD27 = GeographicalCRS("+ellps=clrk66", "NAD27 (Geographical)")
LonLatNAD83 = GeographicalCRS("+ellps=GRS80", "NAD83 (Geographical)")

UPSNorth = ProjectedCRS(proj="+proj=stere +lat_0=90 +lat_ts=90 +lon_0=0 "
                             "+k=0.994 +x_0=2000000 +y_0=2000000 +units=m "
                             "+datum=WGS84 +no_defs",
        spheroid="+ellps=WGS84", name="Universal Polar Stereographic (North)")

UPSSouth = ProjectedCRS(proj="+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=0 "
                             "+k=0.994 +x_0=2000000 +y_0=2000000 +units=m "
                             "+datum=WGS84 +no_defs",
        spheroid="+ellps=WGS84", name="Universal Polar Stereographic (South)")

NSIDCNorth = ProjectedCRS(proj="+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 "
                               "+k=1 +x_0=0 +y_0=0 +units=m +datum=WGS84 "
                               "+no_defs",
        spheroid="+ellps=WGS84", name="NSIDC (North)")

NSIDCSouth = ProjectedCRS(proj="+proj=stere +lat_0=-90 +lat_ts=-70 +lon_0=0 "
                               "+k=1 +x_0=0 +y_0=0 +units=m +datum=WGS84 "
                               "+no_defs",
        spheroid="+ellps=WGS84", name="NSIDC (South)")

LambertEqualArea = ProjectedCRS(proj="+proj=laea +lat_0=0 +lon_0=0 +x_0=0 "
                                     "+y_0=0 +datum=WGS84",
        spheroid="+ellps=WGS84", name="Lambert Equal Area")

GallPetersEqualArea = ProjectedCRS("+proj=cea +lon_0=0 +lat_ts=45 +x_0=0 +y_0=0 "
                                   "+datum=WGS84 +units=m +no_defs",
        spheroid="+ellps=WGS84", name="Gall Peters Equal Area")

WebMercator = ProjectedCRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 "
                           "+lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m "
                           "+nadgrids=@null +wktext +no_defs",
        spheroid="+a=6378137 +b=6378137", name="Web Mercator")

