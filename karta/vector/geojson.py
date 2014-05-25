""" This module includes the `GeoJSON` class and convenience functions for
dealing with GeoJSON data. The `GeoJSON` class uses the builtin json module and
enforces GeoJSON standards. """

import sys
import copy
import json
from collections import namedtuple
import itertools
from numbers import Number
from math import isnan
import traceback


Point = namedtuple('Point', ['coordinates', 'crs'])
MultiPoint = namedtuple('MultiPoint', ['coordinates', 'crs'])
LineString = namedtuple('LineString', ['coordinates', 'crs'])
MultiLineString = namedtuple('MultiLineString', ['coordinates', 'crs'])
Polygon = namedtuple('Polygon', ['coordinates', 'crs'])
MultiPolygon = namedtuple('MultiPolygon', ['coordinates', 'crs'])
GeometryCollection = namedtuple('GeometryCollection', ['geometries'])


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
    """ Class for converting guppy objects to GeoJSON strings. Multipoint-based
    opbjects are written as 'Features'.

    CRS defaults to be named. A linked CRS can be used by passing
    `linkedcrs=True`.

    Notes:
    ------
    Does not handle 'FeatureCollection' types (see printFeatureCollection()
    function).
    """
    def __init__(self, gpobj, linkedcrs=False, **kwargs):
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

        self.gpobj = gpobj
        self.supobj = {}
        if crs is not None:
            self.add_crs(crs, linkedcrs=linkedcrs)
        if self.typestr != 'Point':
            self.supobj['type'] = 'Feature'
            self.add_bbox(bbox)
            self.add_geometry()
        else:
            self.supobj['type'] = 'Point'
            self.add_coordinates()
        return

    def add_crs(self, crs, linkedcrs, **kw):
        """ Tag the JSON object with coordinate reference system metadata using
        a CRS instance, which is expected the have the interface of crs.CRS.

        If True, `linkedcrs` indicates a web-referenced linked CRS. If False,
        the CRS is taked to be named.

        In the case of a linked CRS, the link address and type can be specified
        using the `href` and `linktype` keyword arguments. If not supplied,
        then it is assumed that crs.urn contains the *href* and crs.type
        contains the *type*.

        For more details, see the GeoJSON specification at:
        http://geojson.org/geojson-spec.html#coordinate-reference-system-objects
        """
        if linkedcrs:
            href = kw.get('href', crs.urn)
            linktype = kw.get('linktype', crs.type)
            self.supobj['crs'] = {'type': 'link',
                                  'properties': {'href': href, 'type': linktype}}

        else:
            self.supobj['crs'] = {'type': 'name',
                                  'properties': {'name': crs.urn}}
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
    """ Class for reading general GeoJSON strings and converting to appropriate
    guppy types. Initialize with a file- or StreamIO-like object, and then use
    the pull* methods to extract geometry types. The equivalency table looks
    like:

    -----------------------------------------------------
    | Guppy class       | GeoJSON geometry object       |
    --------------------+--------------------------------
    | Point             | Point                         |
    |                   |                               |
    | Multipoint        | MultiPoint                    |
    |                   |                               |
    | Line              | Linestring, MultiLineString   |
    |                   |                               |
    | Polygon           | Polygon, MultiPolygon         |
    |                   |                               |
    -----------------------------------------------------
     """

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
        self.crs = self.getcrs()
        return

    def _walk(self, dic, geotype):
        """ Find all instances of `key` using recursion. """
        if 'type' in dic and geotype == dic['type']:
            yield dic
        for k in dic:
            if hasattr(dic[k], 'keys'):
                for val in self._walk(dic[k], geotype):
                    yield val
            elif k == 'features':
                for feature in dic[k]:
                    for val in self._walk(feature, geotype):
                        yield val
        return

    def _contains(self, dic, key):
        """ Return whether *key* exists within the JSON hierarchy *dic*. """
        if key in dic:
            return True
        else:
            for val in dic.values():
                if self._contains(val, key):
                    return True
        return False

    def _getproperties(self):
        """ Read feature properties and return a two dictionaries:
            - one representing singleton values ("properties")
            - one representing per-vertex values ("data")
        """
        raise NotImplementedError
        return {}, {}

    def getcrs(self):
        """ Look in the top level for a `crs` field. Return a tuple that is
        either (name, _), (href, type), or (_, _) for a named, linked, or
        undefined CRS, respectively. """
        crs = self.jsondict.get("crs", None)
        if crs is not None:
            if "name" in crs["properties"]:
                res = (crs["properties"]["name"], None)
            else:
                res = (crs["properties"]["href"], crs["properties"]["type"])
        else:
            res = (None, None)
        return res

    def pull_features(self):
        """ Find and return all Feature objects from a FeatureCollection """
        jsonfeats = [obj for obj in self._walk(self.jsondict, 'Feature')]

        featuretype = lambda a: a.get('geometry', {}).get('type', None)
        ispoint = lambda a: featuretype(a) == 'Point'
        ismultipoint = lambda a: featuretype(a) == 'MultiPoint'
        isline = lambda a: featuretype(a) in ('LineString', 'MultiLineString')
        ispolygon = lambda a: featuretype(a) in ('Polygon', 'MultiPolygon')

        # This is an idiotic algorithm - pre-mapping the types isn't saving
        # time because each feature needs to be iterated through anyway
        points = filter(ispoint, jsonfeats)
        multipoints = filter(ismultipoint, jsonfeats)
        lines = filter(isline, jsonfeats)
        polygons = filter(ispolygon, jsonfeats)

        features = []
        for feat in points:
            features.append((self.pull_points(feat)[0], feat['properties'],
                             feat.get('id', None)))
        for feat in multipoints:
            features.append((self.pull_multipoints(feat)[0], feat['properties'],
                             feat.get('id', None)))
        for feat in lines:
            features.append((self.pull_lines(feat)[0], feat['properties'],
                             feat.get('id', None)))
        for feat in polygons:
            features.append((self.pull_polygons(feat)[0], feat['properties'],
                             feat.get('id', None)))

        return features

    def pull_points(self, dic=None):
        """ Return a list of all geometries that can be coerced into a single
        Point. """
        if dic is None:
            dic = self.jsondict
        jsonpoints = self._walk(dic, 'Point')
        points = []
        for point in jsonpoints:
            points.append(Point(point['coordinates'], self.crs))
        return points

    def pull_multipoints(self, dic=None):
        """ Return a list of all geometries that can be coerced into a single
        Multipoint. """
        if dic is None:
            dic = self.jsondict
        jsonmultipoints = self._walk(dic, 'MultiPoint')
        multipoints = []
        for jsonmultipoint in jsonmultipoints:
            multipoints.append(MultiPoint(jsonmultipoint['coordinates'], self.crs))
        return multipoints

    def pull_lines(self, dic=None):
        """ Return a list of all geometries that can be coerced into a Line.
        """
        if dic is None:
            dic = self.jsondict
        jsonlines = self._walk(dic, 'LineString')
        jsonmultilines = self._walk(dic, 'MultiLineString')
        lines = []
        for jsonline in jsonlines:
            lines.append(LineString(jsonline['coordinates'], self.crs))
        for jsonmultiline in jsonmultilines:
            for vertices in jsonmultiline['coordinates']:
                lines.append(MultiLineString(vertices, self.crs))
        return lines

    def pull_polygons(self, dic=None):
        """ Return a list of all geometries that can be coerced into a Polygon.
        """
        if dic is None:
            dic = self.jsondict
        jsonpolygons = self._walk(dic, 'Polygon')
        jsonmultipolygons = self._walk(dic, 'MultiPolygon')
        polygons = []
        for jsonpolygon in jsonpolygons:
            polygons.append(Polygon(jsonpolygon['coordinates'], self.crs))
        for jsonmultipolygon in jsonmultipolygons:
            polygons.append(Polygon(jsonmultipolygon['coordinates'], self.crs))
        return polygons

    def iter_geometries(self):
        """ Return an iterator through all geometries. """
        itergeo = itertools.chain(self.pull_points(),
                                  self.pull_multipoints(),
                                  self.pull_lines(),
                                  self.pull_polygons())
        return itergeo


def list_rec(A):
    """ Recursively convert nested iterables to nested lists """
    if hasattr(A, '__iter__'):
        return [list_rec(el) for el in A]
    else:
        return A

def printGeometryCollection(gpobj_list, **kwargs):
    """ Given an iterable that returns guppy objects, construct a GeoJSON
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


