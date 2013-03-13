""" Provides functions for reading and writing ESRI shapefiles and returning a
guppy object. """

import os
from shapefile import Reader, Writer
import guppy
import traceback

# # Constants for shape types
# NULL = 0
# POINT = 1
# POLYLINE = 3
# POLYGON = 5
# MULTIPOINT = 8
# POINTZ = 11
# POLYLINEZ = 13
# POLYGONZ = 15
# MULTIPOINTZ = 18
# POINTM = 21
# POLYLINEM = 23
# POLYGONM = 25
# MULTIPOINTM = 28
# MULTIPATCH = 31

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

        reader = Reader(shp=files['shp'], shx=files['shx'], dbf=files['dbf'])
        features = []
        for shape in reader.shapes():
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

def write_shapefile(features, stem):
    """ Write a list of features to a shapefile. The features must be of the
    same type. """
    fnms = get_filenames(stem, check=False)

    writer = Writer()
    if not hasattr(features[0], '__iter__'):
        # Must be a Point type
        for feature in features:
            writer.point(*feature.get_vertex())
    else:
        if isinstance(feature, guppy.Line):
            shape_type = 3
        elif isinstance(feature, guppy.Polygon):
            shape_type = 5
        elif isinstance(feature, guppy.Multipoint):
            shape_type = 8
        else:
            raise NotImplementedError("cannot save type "
                                      "{0}".format(type(feature)))

        for feature in features:
            writer.poly([feature.get_vertices()], shape_type)

    try:
        files = open_file_dict(fnms)
        writer.saveShp(files['shp'])
        writer.saveShx(files['shx'])
        writer.saveDbf(files['dbf'])

    finally:
        for f in files.values():
            f.close()

