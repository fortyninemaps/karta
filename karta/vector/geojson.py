""" This module includes the `GeoJSON` class and convenience functions for
dealing with GeoJSON data. The `GeoJSON` class uses the builtin json module and
enforces GeoJSON standards. """

import copy
import json
from collections import namedtuple
from numbers import Number
from math import isnan
from functools import reduce
from ..crs import CRS

Point = namedtuple('Point', ['coordinates', 'crs'])
MultiPoint = namedtuple('MultiPoint', ['coordinates', 'crs'])
LineString = namedtuple('LineString', ['coordinates', 'crs'])
MultiLineString = namedtuple('MultiLineString', ['coordinates', 'crs'])
Polygon = namedtuple('Polygon', ['coordinates', 'crs'])
MultiPolygon = namedtuple('MultiPolygon', ['coordinates', 'crs'])
GeometryCollection = namedtuple('GeometryCollection', ['geometries', 'crs'])
Feature = namedtuple('Feature', ['geometry', 'properties', 'id', 'crs'])
FeatureCollection = namedtuple('FeatureCollection', ['features', 'crs'])

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
        try:
            if o.dtype in ("int8", "int16", "int32", "int64"):
                return int(o)
            elif o.dtype in ("float16", "float32", "float64", "float128"):
                return float(o)
            elif o.dtype in ("complex64", "complex128", "complex256"):
                return complex(o)
            else:
                raise TypeError("not a recognized type")
        except (AttributeError, TypeError):
            return json.JSONEncoder.default(self, o)


class GeoJSONWriter(object):
    """ Class for converting geometry objects to GeoJSON strings. Multipoint-based
    opjects are written as 'Features'.

    Notes:
    ------
    Does not handle 'FeatureCollection' types (see printFeatureCollection()
    function).
    """
    def __init__(self, gpobj, **kwargs):
        type_equiv = {'Point'       : 'Point',
                      'Multipoint'  : 'MultiPoint',
                      'Line'        : 'LineString',
                      'Polygon'     : 'Polygon'}

        if gpobj._geotype in type_equiv:
            self.typestr = type_equiv[gpobj._geotype]
        else:
            raise TypeError('Input object not a recognized geometry')


        self.gpobj = gpobj
        self.supobj = {}

        self.crs = None     # Becomes not None if a non-native CRS is added

        if 'urn' in kwargs:
            self.add_crs(urn=kwargs['urn'])
        elif gpobj.crs is None:
            self.add_crs(crs=DEFAULTCRS)
        else:
            self.add_crs(crs=gpobj.crs)

        if self.typestr != 'Point':
            self.supobj['type'] = 'Feature'
            bbox = kwargs.get('bbox', None)
            self.add_bbox(bbox)
            self.add_geometry()
        else:
            self.supobj['type'] = 'Point'
            self.add_coordinates()
        return

    def add_crs(self, crs=None, urn=None):
        """ Tag the JSON object with coordinate reference system metadata using
        a CRS instance, which is expected the have the interface of crs.CRS.

        In the case of a linked CRS, the link address and type can be specified
        using the `href` and `linktype` keyword arguments. If not supplied,
        then it is assumed that crs.id["urn"] contains the *href* and the
        *type* is "name".

        For more details, see the GeoJSON specification at:
        http://geojson.org/geojson-spec.html#coordinate-reference-system-objects
        """
        if urn is not None:
            self.supobj['crs'] = {'type': 'name',
                                  'properties': {'name': urn}}
        else:
            if hasattr(crs, "jsonhref") and hasattr(crs, "jsontype"):
                self.supobj['crs'] = {'type': 'link',
                                      'properties': {'href': crs.jsonname,
                                                     'type': crs.jsontype}}
            elif hasattr(crs, "jsonname"):
                self.supobj['crs'] = {'type': 'name',
                                      'properties': {'name': crs.jsonname}}
            else:
                self.crs = crs
                self.supobj['crs'] = {'type': 'link',
                                      'properties': {'href': "",
                                                     'type': 'proj4'}}
        return

    def add_bbox(self, bbox=None):
        """ Tag the JSON object with a bounding box. If *bbox* is None
        (default), the bounds will be calculated from the data. Otherwise,
        *bbox* should be in the form (xmin, xmax, ....) for all dimensions.
        """
        if bbox is None:
            if hasattr(self.gpobj, 'get_extent'):
                xmn, xmx, ymn, ymx = self.gpobj.get_extent()
                bbox = ((xmn, xmx), (ymn, ymx))
            else:
                raise AttributeError("Type {0} has no 'get_extent' "
                                     "method".format(type(self.gpobj)))
        self.supobj['bbox'] = bbox
        return

    def add_geometry(self):
        """ Add geometry subobject. This should be called for writing
        Multipoint-derived classes. """
        self.supobj['geometry'] = {'type' : self.typestr}
        self.add_coordinates(self.supobj['geometry'])
        self.add_id()
        self.add_properties()
        return

    def add_coordinates(self, target=None):
        """ Add coordinates. Behaviour depends on the type of the geometry
        object. Target must be None (default) or a dictionary member of the
        JSON object. """
        if target is None:
            target = self.supobj

        if self.gpobj._geotype == 'Polygon':
            target['coordinates'] = [asfixedlist(self.gpobj.get_vertices())]
            if hasattr(self.gpobj, "subs"):
                for poly in self.gpobj.subs:
                    target['coordinates'].append(asfixedlist(poly.get_vertices()))
        elif hasattr(self.gpobj, 'get_vertices'):
            target['coordinates'] = asfixedlist(self.gpobj.get_vertices())
        elif hasattr(self.gpobj, 'get_vertex'):
            target['coordinates'] = asfixedlist(self.gpobj.get_vertex())
        else:
            raise AttributeError('Geometry object has no vertex method')
        return

    def add_properties(self, target=None):
        """ Add the data from the geometry object's data dictionary. """
        if target is None:
            target = self.supobj
        target['properties'] = {}
        if self.gpobj.data is not None:
            if hasattr(self.gpobj.data, 'fields'):
                data = self.gpobj.data
                data = dict((f, data.getfield(f)) for f in data.fields)
            else:
                data = {'point_data': self.gpobj.data}
            for key in data:
                nan2none = lambda a: None if isinstance(a, Number) and isnan(a) else a
                target['properties'][key] = list(map(nan2none, data[key]))
        return

    def add_id(self, target=None):
        """ Add an index field. """
        if target is None:
            target = self.supobj
        target['id'] = list(range(len(self.gpobj.get_vertices())))
        return

    def print_json(self, indent):
        """ Print GeoJSON representation. """
        return json.dumps(self.supobj, indent=indent, cls=NumpyAwareJSONEncoder)

    def write_json(self, fout, indent):
        """ Dump internal dict-object to JSON using the builtin `json` module.
        """
        if self.crs is not None:
            # Indicates that a CRS object should be used to create a linked CRS
            fout_crs = fout.name+".proj4"
            self.supobj["crs"]["properties"]["href"] = fout_crs
            with open(fout_crs, "w") as f:
                f.write(self.crs.get_proj4())
        json.dump(self.supobj, fout, indent=indent, cls=NumpyAwareJSONEncoder)
        return

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
        json_string = serializer.serialize(named_tuple)
    """

    def __init__(self):
        self.enc = NumpyAwareJSONEncoder()
        return

    def serialize(self, geom):
        return self.enc.encode(self.geometry_asdict(geom))

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
        return {"type": name,
                "coordinates": geom.coordinates,
                "crs": self.crsdict(crs=geom.crs)}

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

def printGeometryCollection(gpobj_list, **kwargs):
    """ Given an iterable that returns geometry objects, construct a GeoJSON
    FeatureCollection string.
    """
    geometrylist = []
    for gp in gpobj_list:
        writer = GeoJSONWriter(gp)
        geometrylist.append(copy.copy(writer.supobj))
    geometry_coll = {'type' : 'GeometryCollection',
                    'geometries' : geometrylist}
    kwargs.setdefault("indent", 2)
    return json.dumps(geometry_coll, **kwargs)

def geojson2csv(fin, fout):
    """ Convert a JSON file to a CSV. Only properties (attributes) common to
    all features are retained. The arguments `fin` and `fout` may be either
    filenames or file-like objects.
    """

    if hasattr(fin, 'read'):
        jsdat = json.load(fin)
    else:
        with open(fin, 'r') as f:
            jsdat = json.load(f)

    # Load coordinates and properties
    X = []
    Y = []
    PROP = []
    for feat in jsdat['features']:
        x, y = feat['geometry']['coordinates']
        X.append(x)
        Y.append(y)
        PROP.append(feat['properties'])

    # Find subset of properties that exist in all features
    if len(PROP) > 1:
        s = set(PROP[0].keys())
        common_props = s.intersection(*[set(p.keys()) for p in PROP[1:]])
    else:
        common_props = PROP[0].keys()


    # Output
    outstr = 'x,y' + \
             reduce(lambda a,b: a+','+b, common_props) + \
             '\n'   # Construct header

    outstr += reduce(lambda a, b: a+'\n'+b,                         # Combine lines
                map(lambda x, y, p: str(x) + ',' + str(y) + ',' +   # Construct lines
                    reduce(lambda a, b: str(a) + ',' + str(b),      # Construct properties
                        [p[cp] for cp in common_props]),            # Filter spurious props
                    X, Y, PROP))                                    # Args for line lambda

    try:
        if hasattr(fout, 'read'):
            f = fout
        else:
            f = open(fout, 'w')
        f.write(outstr)

    finally:
        if not hasattr(fout, 'read'):
            f.close()

    return
