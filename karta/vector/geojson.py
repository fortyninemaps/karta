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
GeometryCollection = namedtuple('GeometryCollection', ['geometries'])
Feature = namedtuple('Feature', ['geometry', 'properties', 'id', 'crs'])
FeatureCollection = namedtuple('FeatureCollection', ['features', 'crs'])

class GeoJSONNamedCRS(CRS):
    def __init__(self, name):
        self.name = name
        self.jsonname = name

class GeoJSONLinkedCRS(CRS):
    def __init__(self, href, linktype):
        self.name = href
        self.jsonhref = href
        self.jsontype = linktype

DEFAULTCRS = GeoJSONNamedCRS("urn:ogc:def:crs:epsg::4326")

class ExtendedJSONEncoder(json.JSONEncoder):
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
    opbjects are written as 'Features'.

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
        crs = gpobj._crs
        bbox = kwargs.get('bbox', None)
        urn = kwargs.get('urn', None)

        self.gpobj = gpobj
        self.supobj = {}
        if crs is not None:
            self.add_crs(crs, urn=urn)
        if self.typestr != 'Point':
            self.supobj['type'] = 'Feature'
            self.add_bbox(bbox)
            self.add_geometry()
        else:
            self.supobj['type'] = 'Point'
            self.add_coordinates()
        return

    def add_crs(self, crs, **kw):
        """ Tag the JSON object with coordinate reference system metadata using
        a CRS instance, which is expected the have the interface of crs.CRS.

        In the case of a linked CRS, the link address and type can be specified
        using the `href` and `linktype` keyword arguments. If not supplied,
        then it is assumed that crs.id["urn"] contains the *href* and the
        *type* is "name".

        For more details, see the GeoJSON specification at:
        http://geojson.org/geojson-spec.html#coordinate-reference-system-objects
        """
        if hasattr(crs, "jsonhref") and hasattr(crs, "jsontype"):
            self.supobj['crs'] = {'type': 'link',
                                  'properties': {'href': crs.jsonname,
                                                 'type': crs.jsontype}}
        elif hasattr(crs, "jsonname"):
            self.supobj['crs'] = {'type': 'name',
                                  'properties': {'name': crs.jsonname}}
        else:
            urn = kw.get("urn", geturn(crs))
            self.supobj['crs'] = {'type': 'name',
                                  'properties': {'name': urn}}
        return

    def add_bbox(self, bbox=None):
        """ Tag the JSON object with a bounding box. If *bbox* is None
        (default), the bounds will be calculated from the data. Otherwise,
        *bbox* should be in the form (xmin, xmax, ....) for all dimensions.
        """
        if bbox is None:
            if hasattr(self.gpobj, 'get_extents'):
                bbox = self.gpobj.get_extents()
            else:
                raise AttributeError("Type {0} has no 'get_extents' "
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
            target['coordinates'] = [list_rec(self.gpobj.get_vertices())]
            if hasattr(self.gpobj, "subs"):
                for poly in self.gpobj.subs:
                    target['coordinates'].append(list_rec(poly.get_vertices()))
        elif hasattr(self.gpobj, 'get_vertices'):
            target['coordinates'] = list_rec(self.gpobj.get_vertices())
        elif hasattr(self.gpobj, 'get_vertex'):
            target['coordinates'] = list_rec(self.gpobj.get_vertex())
        else:
            raise AttributeError('Geometry object has no vertex method')
        return

    def add_properties(self, target=None):
        """ Add the data from the geometry object's data dictionary. """
        if target is None:
            target = self.supobj
        target['properties'] = {}
        if (False in (a is None for a in self.gpobj.data)):
            if hasattr(self.gpobj.data, 'keys'):
                data = self.gpobj.data
            else:
                data = {'point_data': self.gpobj.data}
            for key in data:
                nan2none = lambda a: None if isinstance(a, Number) and isnan(a) else a
                target['properties'][key] = list(map(nan2none, data.getfield(key)))
        return

    def add_id(self, target=None):
        """ Add an index field. """
        if target is None:
            target = self.supobj
        target['id'] = list(range(len(self.gpobj.get_vertices())))
        return

    def print_json(self):
        """ Print GeoJSON representation. """
        return json.dumps(self.supobj, indent=2, cls=ExtendedJSONEncoder)

    def write_json(self, fout):
        """ Dump internal dict-object to JSON using the builtin `json` module.
        """
        json.dump(self.supobj, fout, indent=2, cls=ExtendedJSONEncoder)
        return

class GeoJSONReader(object):

    def __init__(self, finput):
        """ Create a reader-object for a GeoJSON-containing file or StreamIO
        object. Use as::

            with open(`fnm`, 'r') as f:
                reader = GeoJSONReader(f)
        """
        if not hasattr(finput, 'read'):
            with open(finput) as f:
                self.jsondict = json.load(f)
        else:
            self.jsondict = json.load(finput)
        return

    @staticmethod
    def parsePoint(d, defaultcrs):
        crs = d.get("crs", defaultcrs)
        return Point(d["coordinates"], crs)

    @staticmethod
    def parseMultiPoint(d, defaultcrs):
        crs = d.get("crs", defaultcrs)
        return MultiPoint(d["coordinates"], crs)

    @staticmethod
    def parseLineString(d, defaultcrs):
        crs = d.get("crs", defaultcrs)
        return LineString(d["coordinates"], crs)

    @staticmethod
    def parseMultiLineString(d, defaultcrs):
        crs = d.get("crs", defaultcrs)
        return MultiLineString(d["coordinates"], crs)

    @staticmethod
    def parsePolygon(d, defaultcrs):
        crs = d.get("crs", defaultcrs)
        return Polygon(d["coordinates"], crs)

    @staticmethod
    def parseMultiPolygon(d, defaultcrs):
        crs = d.get("crs", defaultcrs)
        return MultiPolygon(d["coordinates"], crs)

    def parseGeometry(self, o, defaultcrs):
        crs = o.get("crs", defaultcrs)
        t = o["type"]
        if t == "Point":
            return self.parsePoint(o, crs)
        elif t == "MultiPoint":
            return self.parseMultiPoint(o, crs)
        elif t == "LineString":
            return self.parseLineString(o, crs)
        elif t == "MultiLineString":
            return self.parseMultiLineString(o, crs)
        elif t == "Polygon":
            return self.parsePolygon(o, crs)
        elif t == "MultiPolygon":
            return self.parseMultiPolygon(o, crs)
        else:
            raise TypeError("Unrecognized type {0}".format(t))

    def parseGeometryCollection(self, o, defaultcrs):
        crs = o.get("crs", defaultcrs)
        geoms = [self.parseGeometry(g) for g in o["geometries"]]
        return GeometryCollection(geoms, crs)

    def parseFeature(self, o, defaultcrs):
        crs = o.get("crs", defaultcrs)
        geom = self.parseGeometry(o["geometry"], crs)
        prop = self.parseProperties(o["properties"], len(geom.coordinates))
        fid = o.get("id", None)
        return Feature(geom, prop, fid, crs)

    def parseFeatureCollection(self, o, defaultcrs):
        crs = o.get("crs", defaultcrs)
        features = [self.parseFeature(f, crs) for f in o["features"]]
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

    def items(self, o=None):
        if o is None:
            o = self.jsondict
        crs = o.get("crs", DEFAULTCRS)
        items = []
        if o["type"] == "GeometryCollection":
            items.append(self.parseGeometryCollection(o, crs))
        elif o["type"] == "FeatureCollection":
            items.append(self.parseFeatureCollection(o, crs))
        elif o["type"] == "Feature":
            items.append(self.parseFeature(o, crs))
        else:
            items.append(self.parseGeometry(o, crs))
        return items

def geturn(crs):
    """ Return a urn string for *crs* """
    print("Warning: CRS URN lookup not implemented")
    print(crs)
    print("returning CRS name instead")
    return crs.name


def list_rec(A):
    """ Recursively convert nested iterables to nested lists """
    if hasattr(A, '__iter__'):
        return [list_rec(el) for el in A]
    else:
        return A

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


