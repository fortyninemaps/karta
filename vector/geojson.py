""" This module includes the `GeoJSON` class and convenience functions for
dealing with GeoJSON data. The `GeoJSON` class uses the builtin json module and
enforces GeoJSON standards. """

import guppy
import json
import traceback

class GeoJSONWriter:
    """ Object for converting guppy objects to GeoJSON strings.
    Subclasses must implement:
        - a data writer
        - calc_bbox()
    as well as add an attribute 'strtype'

    DOES NOT HANDLE FEATURE COLLECTIONS
    """
    supobj = {}

    def __init__(self, gpobj, fnm, crs_urn=None, bbox=None):
        type_equiv = {guppy.Point:'Point',
                      guppy.Multipoint:'MultiPoint',
                      guppy.Line:'LineString',
                      guppy.Polygon:'Polygon'}

        if type(gpobj) in type_equiv:
            self.typestr = type_equiv[type(gpobj)]
        else:
            raise NotImplementedError('input object not a recognizes guppy '
                                      'type')
        self.gpobj = gpobj
        if crs_urn is not None:
            self.add_named_crs(crs_urn)
        if self.typestr != 'Point':
            self.supobj['type'] = 'Feature'
            self.add_bbox(bbox)
            self.add_geometry()
        else:
            self.supobj['type'] = 'Point'
            self.add_coordinates()
        self.write_json(fnm)
        return

    def add_named_crs(self, crs_urn):
        """ Tag the JSON object with coordinate reference system metadata using
        a named CRS. The crs_urn should be a valid OGC CRS URN identifier. """
        # Validate the crs_urn?
        supobj['crs'] = {'type'         :   'name',
                         'properties'   :   {'name' :   crs_urn}}
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
        if hasattr(self.gpobj.data, 'keys') or \
                (False in (a is None for a in self.gpobj.data)):
            self.add_properties()
        return

    def add_coordinates(self, target=None):
        """ Add coordinates. Behaviour depends on the type of the geometry
        object. Target must be None (default) or a dictionary member of the
        JSON object. """
        if target is None:
            target = self.supobj

        if hasattr(self.gpobj, 'get_vertices'):
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
        if hasattr(self.gpobj.data, 'keys'):
            data = self.gpobj.data
        else:
            data = {'point_data': self.gpobj.data}
        target['properties'] = {}
        for key in data.keys():
            target['properties'][key] = data[key]
        return

    def add_id(self, target=None):
        """ Add an index field. """
        if target is None:
            target = self.supobj
        target['id'] = range(len(self.gpobj.get_vertices()))
        return

    def write_json(self, fnm):
        """ Dump internal dict-object to JSON using the builtin `json` module. """
        with open(fnm, 'w') as f:
            json.dump(self.supobj, f, indent=2)
        return

#class GeoJSONReader:


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
    outstr = 'x,y' + reduce(lambda a,b: a+','+b, common_props) + '\n'   # Construct header

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


