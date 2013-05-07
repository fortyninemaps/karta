""" This module includes the `GeoJSON` class and convenience functions for
dealing with GeoJSON data. The `GeoJSON` class uses the builtin json module and
enforces GeoJSON standards. """

import sys
import copy
import guppy
import json
import traceback

class GeoJSONWriter(object):
    """ Class for converting guppy objects to GeoJSON strings. Multipoint-based
    opbjects are written as 'Features'.

    Notes:
    ------
    Does not handle 'FeatureCollection' types (see printFeatureCollection()
    function).
    """
    supobj = {}

    def __init__(self, gpobj, **kwargs):
        type_equiv = {guppy.Point       :   'Point',
                      guppy.Multipoint  :   'MultiPoint',
                      guppy.Line        :   'LineString',
                      guppy.Polygon     :   'Polygon'}

        if type(gpobj) in type_equiv.keys():
            self.typestr = type_equiv[type(gpobj)]
        else:
            raise NotImplementedError('input object not a recognizes guppy '
                                      'type')
        crs = kwargs.get('crs', None)
        crs_fmt = kwargs.get('crs_fmt', None)
        bbox = kwargs.get('bbox', None)

        self.gpobj = gpobj
        if crs is not None:
            self.add_named_crs(crs, crs_fmt)
        if self.typestr != 'Point':
            self.supobj['type'] = 'Feature'
            self.add_bbox(bbox)
            self.add_geometry()
        else:
            self.supobj['type'] = 'Point'
            self.add_coordinates()
        return

    def add_named_crs(self, crs, fmt='epsg'):
        """ Tag the JSON object with coordinate reference system metadata using
        a named CRS. `fmt` indicates the format of the named CRS argument, and
        may be one of ('epsg', 'ogc_crs_urn').
        """
        if fmt == 'epsg':
            href = ('http://spatialreference.org/ref/epsg/{0}/'
                    'proj4/'.format(crs))
            islinked = True
        elif fmt == 'ogc_crs_urn':
            crs_str = crs
            islinked = False
        else:
            raise NotImplementedError('CRS fmt {0} not recognized'.format(fmt))

        if islinked:
            self.supobj['crs'] = {'type': 'link',
                                  'properties': {'href': href, 'type': 'proj4'}}
        else:
            self.supobj['crs'] = {'type': 'name',
                                  'properties': {'name': crs_str}}
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
        self.supobj['bbox'] = {'bbox'    :   bbox}
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

        if isinstance(self.gpobj, guppy.Polygon):
            target['coordinates'] = [self.gpobj.get_vertices()]
            if hasattr(self.gpobj, "subs"):
                for poly in self.gpobj.subs:
                    target['coordinates'].append(poly.get_vertices())
        elif hasattr(self.gpobj, 'get_vertices'):
            target['coordinates'] = self.gpobj.get_vertices()
        elif hasattr(self.gpobj, 'get_vertex'):
            target['coordinates'] = self.gpobj.get_vertex()
        else:
            raise guppy.GGeoError('Geometry object has no vertex method')
        return

    def add_properties(self, target=None):
        """ Add the data from the geometry object's data attribute. """
        if target is None:
            target = self.supobj
        target['properties'] = {}
        if (False in (a is None for a in self.gpobj.data)):
            if hasattr(self.gpobj.data, 'keys'):
                data = self.gpobj.data
            else:
                data = {'point_data': self.gpobj.data}
            for key in data.keys():
                target['properties'][key] = list(data[key])
        return

    def add_id(self, target=None):
        """ Add an index field. """
        if target is None:
            target = self.supobj
        target['id'] = range(len(self.gpobj.get_vertices()))
        return

    def print_json(self):
        """ Print GeoJSON representation. """
        return json.dumps(self.supobj, indent=2)

    def write_json(self, fout):
        """ Dump internal dict-object to JSON using the builtin `json` module.
        """
        json.dump(self.supobj, fout, indent=2)
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

    def pull_points(self):
        """ Return a list of all geometries that can be coerced into a single
        Point. """
        jsonpoints = self._walk(self.jsondict, 'Point')
        points = []
        for jsonpoint in jsonpoints:
            points.append(guppy.Point(jsonpoint['coordinates']))
        return points

    def pull_multipoints(self):
        """ Return a list of all geometries that can be coerced into a single
        Multipoint. """
        jsonmultipoints = self._walk(self.jsondict, 'MultiPoint')
        multipoints = []
        for jsonmultipoint in jsonmultipoints:
            multipoints.append(guppy.Multipoint(jsonmultipoint['coordinates']))
        return multipoints

    def pull_lines(self):
        """ Return a list of all geometries that can be coerced into a Line.
        """
        jsonlines = self._walk(self.jsondict, 'LineString')
        jsonmultilines = self._walk(self.jsondict, 'MultiLineString')
        lines = []
        for jsonline in jsonlines:
            lines.append(guppy.Line(jsonline['coordinates']))
        for jsonmultiline in jsonmultilines:
            for vertices in jsonmultiline['coordinates']:
                lines.append(guppy.Line(vertices))
        return lines

    def pull_polygons(self):
        """ Return a list of all geometries that can be coerced into a Polygon.
        """

        def add_polygon(L, parts):
            """ Convert a list of polygon parts to a guppy.Polygon and append
            to a list `L`. """
            if len(parts) > 1:
                raise NotImplementedError('polygon holes')
            L.append(guppy.Polygon(parts[0]))
            return

        jsonpolygons = self._walk(self.jsondict, 'Polygon')
        jsonmultipolygons = self._walk(self.jsondict, 'MultiPolygon')
        polygons = []
        for jsonpolygon in jsonpolygons:
            add_polygon(polygons, jsonpolygon['coordinates'])
        for jsonmultipolygon in jsonmultipolygons:
            for parts in jsonmultipolygon['coordinates']:
                add_polygon(polygons, parts)
        return polygons

    def list_features(self):
        """ List the features present. """
        raise NotImplementedError


def print_FeatureCollection(gpobj_list, **kwargs):
    """ Given an iterable that returns guppy objects, construct a GeoJSON
    FeatureCollection string.
    """
    featurelist = []
    for gp in gpobj_list:
        writer = GeoJSONWriter(gp)
        featurelist.append(copy.copy(writer.supobj))
    feature_coll = {'type' : 'FeatureCollection',
                    'features' : featurelist}
    return json.dumps(feature_coll, indent=2)



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
    try:
        for feat in jsdat['features']:
            x, y = feat['geometry']['coordinates']
            X.append(x)
            Y.append(y)
            PROP.append(feat['properties'])
    except KeyError:
        sys.stderr.write('invalid geojson schema')
        traceback.print_exc()
    except:
        traceback.print_exc()

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


