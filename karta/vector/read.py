""" Convenience reader functions """

import os
import shapefile
import dateutil.parser
from . import guppy
from . import geojson
from . import xyfile
from .. import crs
from .metadata import Metadata

def _parsegeojsoncrs(crstup):
    """ From a tuple representing a GeoJSON (name,None) or (href,type) pair,
    return an appropriate karta crs instance. """
    if crstup[0] is None:       # use default as defined by spec
        return crs.LONLAT_WGS84
    elif crstup[1] is None:     # named CRS
        for c in crs.crslist:
            if c.urn == crstup[0]:
                return c
        return crs.CRS("unknown", "unknown", crstup[0])
    else:                       # linked CRS
        return crs.CRS("unknown", crstup[1], crstup[0])

def read_geojson(f):
    """ Read a GeoJSON object file and return a list of geometries """
    R = geojson.GeoJSONReader(f)
    geoms = R.iter_geometries()
    gplist = []
    for geom in geoms:
        coordsys = _parsegeojsoncrs(geom.crs)
        if isinstance(geom, geojson.Point):
            gplist.append(guppy.Point(geom.coordinates, crs=coordsys))
        elif isinstance(geom, geojson.MultiPoint):
            gplist.append(guppy.Multipoint(geom.coordinates, crs=coordsys))
        elif isinstance(geom, geojson.LineString):
            gplist.append(guppy.Line(geom.coordinates, crs=coordsys))
        elif isinstance(geom, geojson.Polygon):
            gplist.append(guppy.Polygon(geom.coordinates[0],
                                        subs=geom.coordinates[1:],
                                        crs=coordsys))
    return gplist

def read_geojson_features(f):
    """ Read a GeoJSON object file and return a list of features """
    R = geojson.GeoJSONReader(f)
    features = R.pull_features()
    print()
    geoms = []
    coordsys = _parsegeojsoncrs(R.getcrs())
    for (geom, properties, id) in features:

        if isinstance(geom, geojson.Point):
            (p, d) = _geojson_properties2guppy(properties, 1)
            geoms.append(guppy.Point(geom.coordinates, properties=p, data=d,
                                     crs=coordsys))

        elif isinstance(geom, geojson.MultiPoint):
            (p, d) = _geojson_properties2guppy(properties, len(geom.coordinates))
            geoms.append(guppy.Multipoint(geom.coordinates, properties=p, data=d,
                                          crs=coordsys))

        elif isinstance(geom, geojson.LineString):
            (p, d) = _geojson_properties2guppy(properties, len(geom.coordinates))
            geoms.append(guppy.Line(geom.coordinates, properties=p, data=d,
                                    crs=coordsys))

        elif isinstance(geom, geojson.Polygon):
            (p, d) = _geojson_properties2guppy(properties, len(geom.coordinates))
            geoms.append(guppy.Polygon(geom.coordinates[0], properties=p, data=d,
                                       subs=geom.coordinates[1:],
                                       crs=coordsys))
    return geoms

def _geojson_properties2guppy(properties, n):
    """ Takes a dictionary (derived from a GeoJSON properties object) and
    divides it into singleton properties and *n*-degree data. """
    props = {}
    data = {}
    for (key, value) in properties.items():
        if isinstance(value, list) or isinstance(value, tuple):
            if len(value) == n:
                data[key] = value
            else:
                raise ValueError("properties must be singleton or per-vertex")
        else:
            props[key] = value
    return props, data


def read_xyfile(f, delimiter='', header_rows=0, astype=guppy.Multipoint, coordrank=2):
    """ Read an ASCII delimited table and return a guppy object given by *astype*.
    """
    dat = xyfile.load_xy(f, delimiter=delimiter, header_rows=header_rows)
    ncols = dat.shape[1]
    if ncols >= coordrank:
        coords = dat[:,:coordrank]
        if ncols > coordrank:
            data = dat[:,coordrank:]
        else:
            data = None
        return astype(coords, data=data)
    else:
        raise IOError('data table has insufficient number of columns')

### Shapefile functions ###

def get_filenames(stem, check=False):
    """ Given a filename basename, return the associated shapefile paths. If
    `check` is True, ensure that the files exist."""
    shp = stem + '.shp'
    shx = stem + '.shx'
    dbf = stem + '.dbf'
    if check:
        for fnm in (shp, shx, dbf):
            if not os.path.isfile(fnm):
                raise Exception('missing {0}'.format(fnm))
    return {'shp':shp, 'shx':shx, 'dbf':dbf}

def open_file_dict(fdict):
    """ Open each file in a dictionary of filenames and return a matching
    dictionary of the file objects. """
    files = {}
    for ext in fdict.keys():
        files[ext] = open(fdict[ext], 'rb')
    return files

dBase_type_dict = {"I": int,
                   "O": float,
                   "C": str,
                   "@": dateutil.parser.parse,
                   "L": bool}

def recordsasdata(reader):
    """ Interpret shapefile records as a Metadata object """
    d = {}
    idfunc = lambda a: a
    records = [rec for rec in reader.records()]
    for (i,k) in enumerate(reader.fields[1:]):
        f = dBase_type_dict.get(k[1], idfunc)
        d[k[0]] = [f(rec[i]) for rec in records]
    return Metadata(d)

def recordsasproperties(reader):
    """ Interpret shapefile records as a list of properties dictionaries """
    proplist = []
    keys = reader.fields
    idfunc = lambda a: a
    for (i,rec) in enumerate(reader.records()):
        properties = {}
        for (k,v) in zip(keys, rec):
            f = dBase_type_dict.get(k[1], idfunc)
            properties[k[0]] = f(v)
        proplist.append(properties)
    return proplist

def read_shapefile(stem):
    """ Read a shapefile given `stem`, which is the name without an extension.
    """
    fnms = get_filenames(stem, check=True)

    try:
        files = open_file_dict(fnms)
        reader = shapefile.Reader(shp=files['shp'], shx=files['shx'],
                                  dbf=files['dbf'])

        if reader.shapeType == 1:       # Points
            verts = [shp.points[0] for shp in reader.shapes()]
            d = recordsasdata(reader)
            geoms = [guppy.Multipoint(verts, data=d)]

        elif reader.shapeType == 3:     # Lines
            plist = recordsasproperties(reader)
            geoms = []
            for (shp,prop) in zip(reader.shapes(), plist):
                geoms.append(guppy.Line(shp.points, properties=prop))

        elif reader.shapeType == 5:     # Polygon
            plist = recordsasproperties(reader)
            geoms = []
            for (shp,prop) in zip(reader.shapes(), plist):
                geoms.append(guppy.Polygon(shp.points, properties=prop))

        else:
            raise NotImplementedError("Shapefile shape type {0} not "
                "implemented".format(reader.shapeType))

    finally:
        for f in files.values():
            f.close()

    return geoms

