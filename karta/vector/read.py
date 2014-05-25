""" Convenience reader functions """

import os
import shapefile
from . import guppy
from . import geojson
from . import xyfile
from .. import crs

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


def shape2point(shape):
    """ Convert a shapefile._Shape `shape` to a guppy.Point. """
    return guppy.Point(*shape.points)

def shape2line(shape):
    """ Convert a shapefile._Shape `shape` to a guppy.Line. """
    verts = shape.points
    return guppy.Line(verts)

def shape2poly(shape):
    """ Converts a shapefile._Shape `shape` to a guppy.Polygon. """
    verts = shape.points
    return guppy.Polygon(verts)

def shape2multipoint(shape):
    """ Converts a shapefile._Shape `shape` to a guppy.Polygon. """
    verts = shape.points
    return guppy.Multipoint(verts)

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
        files[ext] = open(fdict[ext], 'r')
    return files

def read_shapefile(stem):
    """ Read a shapefile given `stem`, which is the name without an extension.
    """
    fnms = get_filenames(stem, check=True)

    try:
        files = open_file_dict(fnms)

        reader = shapefile.Reader(shp=files['shp'], shx=files['shx'],
                                  dbf=files['dbf'])
        features = []
        for shape in reader.shapes():
            if len(shape.points) > 0:
                if shape.shapeType == 1:
                    features.append(shape2point(shape))
                elif shape.shapeType == 3:
                    features.append(shape2line(shape))
                elif shape.shapeType == 5:
                    features.append(shape2poly(shape))
                elif shape.shapeType == 8:
                    features.append(shape2multipoint(shape))
                else:
                    raise NotImplementedError("cannot read shape type "
                                              "{0}".format(shape.shapeType))

    finally:
        for f in files.values():
            f.close()

    return features

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

