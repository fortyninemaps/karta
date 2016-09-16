"""
GeoJSON drivers for Karta

Defines named tuples representing GeoJSON entities.

Overview
--------

`GeoJSONReader` converts GeoJSON strings to named tuples.

`GeoJSONSerializer` converts named tuples to GeoJSON strings.

`as_named_tuple` function converts karta.geometry classes to equivalent named
tuples.
"""

import json
from collections import namedtuple
from functools import reduce
from .utilities import _reproject, _reproject_nested, _as_nested_lists
from ..crs import CRS, LonLatWGS84

Point = namedtuple('Point', ['coordinates', 'crs'])
MultiPoint = namedtuple('MultiPoint', ['coordinates', 'crs'])
LineString = namedtuple('LineString', ['coordinates', 'crs'])
MultiLineString = namedtuple('MultiLineString', ['coordinates', 'crs'])
Polygon = namedtuple('Polygon', ['coordinates', 'crs'])
MultiPolygon = namedtuple('MultiPolygon', ['coordinates', 'crs'])
GeometryCollection = namedtuple('GeometryCollection', ['geometries', 'crs'])
Feature = namedtuple('Feature', ['geometry', 'properties', 'id', 'crs'])
FeatureCollection = namedtuple('FeatureCollection', ['features', 'crs'])

def as_named_tuple(*geoms, **kwargs):
    """ Convert one or more Geometry instances to GeoJSON-structured named tuples.

    Parameters
    ----------
    *geoms : subtypes of karta.Geometry
        karta vector geometries to convert
    urn : str, optional
        URN defining a specific CRS to use

    Returns
    -------
    Either a Feature or a FeatureCollection

    Raises
    ------
    TypeError
        if one or more of geoms has an unrecognized `_geotype` attribute
    """
    if kwargs.get("urn", None) is not None:
        crs = crs_from_urn(kwargs["urn"])
    else:
        crs = crs_from_karta(geoms[0].crs)

    features = []
    for geom in geoms:
        if geom._geotype == "Point":
            g = Point(geom.vertex, crs)
        elif geom._geotype == "Line":
            g = LineString(_as_nested_lists(geom.vertices), crs)
        elif geom._geotype == "Polygon":
            verts = [_as_nested_lists(geom.vertices_ring)]
            for sub in geom.subs:
                verts.append(_as_nested_lists(sub.vertices))
            g = Polygon(verts, crs)
        elif geom._geotype == "Multipoint":
            g = MultiPoint(_as_nested_lists(geom.vertices), crs)
        elif geom._geotype == "Multiline":
            g = MultiLineString(_as_nested_lists(geom.vertices), crs)
        elif geom._geotype == "Multipolygon":
            g = MultiPolygon(_as_nested_lists(geom.vertices_ring), crs)
        else:
            raise TypeError("unhandled type: {0}".format(type(geom)))

        data = {}
        if hasattr(geom, "data"):
            for field in geom.data.fields:
                data[field] = geom.d[field]

        properties = geom.properties
        properties.update(data)
        features.append(Feature(g, properties, 0, crs))

    if len(features) == 1:
        return features[0]
    else:
        return FeatureCollection(features, crs)

def crs_from_urn(urn):
    return {'type': 'name', 'properties': {'name': urn}}

def crs_from_karta(crs):
    if hasattr(crs, "jsonhref") and hasattr(crs, "jsontype"):
        return {'type': 'link', 'properties': {'href': crs.jsonname, 'type': crs.jsontype}}
    elif hasattr(crs, "jsonname"):
        return {'type': 'name', 'properties': {'name': crs.jsonname}}
    else:
        # NOTE: make sure CRS gets written to file by caller
        return {'type': 'link', 'properties': {'href': "", 'type': 'proj4'}}

class GeoJSONNamedCRS(CRS):
    def __init__(self, name):
        self.name = name
        self.jsonname = name

    def __eq__(self, other):
        return self.jsonname == other.jsonname

class GeoJSONLinkedCRS(CRS):
    def __init__(self, href, linktype):
        self.name = href
        self.jsonhref = href
        self.jsontype = linktype

    def __eq__(self, other):
        return self.jsonhref == other.jsonhref

DEFAULTCRS = {"type": "name",
              "properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84"}}

class NumpyAwareJSONEncoder(json.JSONEncoder):
    """ Numpy-specific numbers prior to 1.9 don't inherit from Python numeric
    ABCs. This class is a hack to coerce numpy values into Python types for
    JSON serialization. """

    def default(self, o):
        if not hasattr(o, "dtype") or (hasattr(o, "__len__") and (len(o) != 1)):
            return json.JSONEncoder.default(self, o)
        elif o.dtype in ("int8", "int16", "int32", "int64"):
            return int(o)
        elif o.dtype in ("float16", "float32", "float64", "float128"):
            return float(o)
        elif o.dtype in ("complex64", "complex128", "complex256"):
            return complex(o)
        else:
            raise TypeError("not a recognized type")

class GeoJSONReader(object):

    def __init__(self, finput, defaultcrs=DEFAULTCRS):
        """ Create a reader-object for a GeoJSON-containing file or StreamIO
        object. Use as::

            with open(`fnm`, 'r') as f:
                reader = GeoJSONReader(f)
        """
        if hasattr(finput, 'read'):
            self.jsondict = json.load(finput)
        elif isinstance(finput, str):
            try:
                self.jsondict = json.loads(finput)
            except ValueError:
                with open(finput) as f:
                    self.jsondict = json.load(f)
        else:
            raise TypeError("*finput* must be a file object, a JSON string, or "
                            "a filename string")

        self.defaultcrs = defaultcrs
        return

    def parsePoint(self, d):
        crs = d.get("crs", self.defaultcrs)
        return Point(d["coordinates"], crs)

    def parseMultiPoint(self, d):
        crs = d.get("crs", self.defaultcrs)
        return MultiPoint(d["coordinates"], crs)

    def parseLineString(self, d):
        crs = d.get("crs", self.defaultcrs)
        return LineString(d["coordinates"], crs)

    def parseMultiLineString(self, d):
        crs = d.get("crs", self.defaultcrs)
        return MultiLineString(d["coordinates"], crs)

    def parsePolygon(self, d):
        crs = d.get("crs", self.defaultcrs)
        return Polygon(d["coordinates"], crs)

    def parseMultiPolygon(self, d):
        crs = d.get("crs", self.defaultcrs)
        return MultiPolygon(d["coordinates"], crs)

    def parseGeometry(self, o):
        t = o["type"]
        if t == "Point":
            return self.parsePoint(o)
        elif t == "MultiPoint":
            return self.parseMultiPoint(o)
        elif t == "LineString":
            return self.parseLineString(o)
        elif t == "MultiLineString":
            return self.parseMultiLineString(o)
        elif t == "Polygon":
            return self.parsePolygon(o)
        elif t == "MultiPolygon":
            return self.parseMultiPolygon(o)
        else:
            raise TypeError("Unrecognized type {0}".format(t))

    def parseGeometryCollection(self, o):
        crs = o.get("crs", self.defaultcrs)
        geoms = [self.parse(g) for g in o["geometries"]]
        return GeometryCollection(geoms, crs)

    def parseFeature(self, o):
        crs = o.get("crs", self.defaultcrs)
        geom = self.parse(o["geometry"])
        if isinstance(geom, GeometryCollection):
            n = 1
        else:
            n = len(geom.coordinates)
        prop = self.parseProperties(o["properties"], n)
        fid = o.get("id", None)
        return Feature(geom, prop, fid, crs)

    def parseFeatureCollection(self, o):
        crs = o.get("crs", self.defaultcrs)
        features = [self.parseFeature(f) for f in o["features"]]
        return FeatureCollection(features, crs)

    @staticmethod
    def parseProperties(prop, geomlen):
        d = {"scalar":{},
             "vector":{}}
        for key, value in prop.items():
            if geomlen > 1 and hasattr(value, "__iter__") and len(value) == geomlen:
                d["vector"][key] = value
            else:
                d["scalar"][key] = value
        return d

    def parse(self, d=None):
        if d is None:
            d = self.jsondict
        if d["type"] == "GeometryCollection":
            res = self.parseGeometryCollection(d)
        elif d["type"] == "FeatureCollection":
            res = self.parseFeatureCollection(d)
        elif d["type"] == "Feature":
            res = self.parseFeature(d)
        else:
            res = self.parseGeometry(d)
        return res

class GeoJSONSerializer(object):
    """ Class for converting GeoJSON named tuples to GeoJSON.

    Usage:

        serializer = GeoJSONSerializer()
        json_string = serializer(named_tuple)
    """

    def __init__(self):
        # Previously used encoder directly, but for now calling json.dumps
        #self.enc = NumpyAwareJSONEncoder()
        return

    def __call__(self, geom, indent=None):
        #return self.enc.encode(self.geometry_asdict(geom), indent=indent)
        return json.dumps(self.geometry_asdict(geom), indent=indent,
                          cls=NumpyAwareJSONEncoder)

    def geometry_asdict(self, geom):

        if hasattr(geom, "properties"):
            return self.feature_asdict(geom)

        elif hasattr(geom, "geometries"):
            return self.geometry_collection_asdict(geom)

        elif hasattr(geom, "features"):
            return self.feature_collection_asdict(geom)

        elif isinstance(geom, MultiPoint):
            return self._geometry_asdict(geom, "MultiPoint")

        elif isinstance(geom, Point):
            return self._geometry_asdict(geom, "Point")

        elif isinstance(geom, MultiLineString):
            return self._geometry_asdict(geom, "MultiLineString")

        elif isinstance(geom, LineString):
            return self._geometry_asdict(geom, "LineString")

        elif isinstance(geom, MultiPolygon):
            return self._geometry_asdict(geom, "MultiPolygon")

        elif isinstance(geom, Polygon):
            return self._geometry_asdict(geom, "Polygon")
        else:
            raise TypeError("cannot serialize type '{0}'".format(type(geom)))

    @staticmethod
    def crsdict(crs=None, urn=None, href="", linktype="proj4"):
        """ Return a dictionary that can be serialized as a GeoJSON coordinate
        system mapping using a *urn* or a *crs* instance.

        In the case of a linked CRS, the link address and type can be specified
        using the `href` and `linktype` keyword arguments.

        For more details, see the GeoJSON specification at:
        http://geojson.org/geojson-spec.html#coordinate-reference-system-objects
        """
        if urn is not None:
            return {'type': 'name',
                    'properties': {'name': urn}}
        elif crs is not None:
            if hasattr(crs, "jsonhref") and hasattr(crs, "jsontype"):
                d = {'type': 'link',
                     'properties': {'href': crs.jsonname,
                                    'type': crs.jsontype}}
            elif hasattr(crs, "jsonname"):
                d = {'type': 'name',
                     'properties': {'name': crs.jsonname}}
            else:
                d = {'type': 'link',
                     'properties': {'href': href,
                                    'type': linktype}}
            return d
        else:
            raise ValueError("either a urn must be specified or a CRS object "
                             "provided")
        return

    def _geometry_asdict(self, geom, name):
        if not isinstance(geom.crs, dict):
            crs = self.crsdict(geom.crs)
        else:
            crs = geom.crs
        return {"type": name,
                "coordinates": geom.coordinates,
                "crs": crs}

    def feature_asdict(self, feature):
        return {"type": "Feature",
                "geometry": self.geometry_asdict(feature.geometry),
                "properties": feature.properties}

    def geometry_collection_asdict(self, geometry_collection):
        return {"type": "GeometryCollection",
                "geometries": [self.geometry_asdict(g)
                               for g in geometry_collection.geometries]}

    def feature_collection_asdict(self, feature_collection):
        return {"type": "FeatureCollection",
                "features": [self.feature_asdict(f)
                             for f in feature_collection.features]}

def asfixedlist(A):
    """ Recursively convert nested iterables or coordinates to nested lists at
    fixed precision. """
    if hasattr(A, '__iter__'):
        return [asfixedlist(el) for el in A]
    else:
        return round(A, 6)

class GeoJSONOutMixin(object):
    """ Mixin class to be added to geometry objects, adding geojson
    functionality.
    """

    _geojson_serializer = GeoJSONSerializer()

    def as_geojson(self, indent=None, urn=None, force_wgs84=True):
        """ Output representation of internal data as a GeoJSON string.

        Parameters
        ----------
        indent : int, optional
            indentation of generated GeoJSON (default 2)
        force_wgs84 : bool, optional
            Forces output to use geographical coordinates with the WGS84 datum,
            as recommended by the GeoJSON draft specification
            (https://datatracker.ietf.org/doc/draft-ietf-geojson/).
            If *urn* is not set, "urn:ogc:def:crs:OGC:1.3:CRS84" is used.
            (default True)
        urn : str, optional
            overrides GeoJSON CRS with provided URN string
        """
        if force_wgs84 and (urn is None):
            urn = "urn:ogc:def:crs:OGC:1.3:CRS84"
        if force_wgs84 and (self.crs != LonLatWGS84):
            kw = dict(properties=self.properties, crs=LonLatWGS84)
            if hasattr(self, "data"):
                kw["data"] = self.data
            if hasattr(self, "vertices"):
                vertices = _reproject_nested(self.vertices, self.crs, LonLatWGS84)
            else:
                vertices = _reproject(self.vertex, self.crs, LonLatWGS84)
            geo = type(self)(vertices, **kw)
        else:
            geo = self

        return self._geojson_serializer(as_named_tuple(geo, urn=urn),
                                        indent=indent)

    def to_geojson(self, f, indent=None, urn=None, force_wgs84=True):
        """ Write data as a GeoJSON string to a file-like object `f`.

        Parameters
        ----------
        f : str or file-like object
            file to receive GeoJSON string
        indent : int, optional
            indentation of generated GeoJSON (default None)
        force_wgs84 : bool, optional
            Forces output to use geographical coordinates with the WGS84 datum,
            as recommended by the GeoJSON draft specification
            (https://datatracker.ietf.org/doc/draft-ietf-geojson/).
            If *urn* is not set, "urn:ogc:def:crs:OGC:1.3:CRS84" is used.
            (default True)
        urn : str, optional
            overrides GeoJSON CRS with provided URN string
        """
        try:
            if not hasattr(f, "write"):
                fobj = open(f, "w")
            else:
                fobj = f
            fobj.write(self.as_geojson(indent=indent, urn=urn,
                                       force_wgs84=force_wgs84))
        finally:
            if not hasattr(f, "write"):
                fobj.close()
        return

