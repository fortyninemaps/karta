""" Convenience reader functions """

import os
import shapefile
from . import guppy
from . import geojson
from . import xyfile

def read_geojson(f):
    """ Read a GeoJSON object file and return a list of geometries """
    R = geojson.GeoJSONReader(f)
    geoms = R.iter_geometries()
    gplist = []
    for geom in geoms:
        if isinstance(geom, geojson.Point):
            gplist.append(guppy.Point(geom.coordinates))
        elif isinstance(geom, geojson.MultiPoint):
            gplist.append(guppy.Multipoint(geom.coordinates))
        elif isinstance(geom, geojson.LineString):
            gplist.append(guppy.Line(geom.coordinates))
        elif isinstance(geom, geojson.Polygon):
            gplist.append(guppy.Polygon(geom.coordinates[0],
                                        subs=geom.coordinates[1:]))
    return gplist

def read_geojson_features(f):
    """ Read a GeoJSON object file and return a list of features """
    R = geojson.GeoJSONReader(f)
    feats = R.pull_features()
    gplist = []
    raise NotImplementedError
    for feat in feats:

        # [...]
        pass

    return gplist


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

