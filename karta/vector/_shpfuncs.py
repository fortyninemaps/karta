""" Provides functions for reading and writing ESRI shapefiles and returning a
guppy object. """

import os
import datetime
import numbers
import shapefile
from shapefile import Reader, Writer
#from . import guppy

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

def open_file_dict(fdict, mode='rb'):
    """ Open each file in a dictionary of filenames and return a matching
    dictionary of the file objects. """
    files = {}
    for ext in fdict.keys():
        files[ext] = open(fdict[ext], mode)
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

def property_field_type(value):
    """ Determine the appropriate dBase field type for *value* """
    if isinstance(value, numbers.Number):
        if isinstance(value, numbers.Integral):
            desc = "I"
        elif isinstance(value, numbers.Real):
            desc = "O"
        else:
            raise TypeError("cannot choose the correct dBase type for {0}\n".format(type(value)))
    elif isinstance(value, str):
        desc = "C"
    elif isinstance(value, datetime.datetime):
        desc = "@"
    elif isinstance(value, datetime.date):
        desc = "D"
    elif isinstance(value, bool):
        desc = "L"
    else:
        raise TypeError("cannot choose the correct dBase type for {0}\n".format(type(value)))
    return desc

def addfields(writer, properties):
    if len(properties) == 0:
        writer.field("PROPERTIES", "C", "40")
        writer.record("_karta.Line")
    else:
        values = []
        for key in properties:
            value = properties[key]
            typ = property_field_type(value)
            length = "100"
            writer.field(key.upper(), typ, length)
            value.append(value)
        writer.record(*values)
    return

def write_multipoint(mp, fstem):
    w = shapefile.Writer(shapeType=shapefile.POINT)
    for pt in mp:
        w.point(*pt.vertex)
    addfields(w, mp.properties)
    w.save(fstem)
    return

def write_line(line, fstem):
    w = shapefile.Writer(shapeType=shapefile.POLYLINE)
    w.poly(shapeType=shapefile.POLYLINE, parts=[line.vertices])
    addfields(w, line.properties)
    w.save(fstem)
    return

def write_poly2(poly, fstem):
    w = shapefile.Writer(shapeType=shapefile.POLYGON)
    w.poly(shapeType=shapefile.POLYGON, parts=[poly.vertices])
    addfields(w, poly.properties)
    w.save(fstem)
    return

def list_parts(feature):
    """ Return a list of polygon parts. """
    return [feature.vertices()].extend(
            [list_parts(p) for p in feature.subs])

def write_shapefile(features, stem):
    """ Write a list of features to a shapefile. The features must be of the
    same type. """
    fnms = get_filenames(stem, check=False)

    writer = Writer()
    if not hasattr(features[0], '__iter__'):
        # Must be a Point type
        for feature in features:
            writer.shapeType = POINT
            writer.point(*feature.get_vertex())
    else:
        if isinstance(features[0], guppy.Line):
            writer.shapeType = shapefile.POLYLINE
        elif isinstance(features[0], guppy.Polygon):
            writer.shapeType = shapefile.POLYGON
        elif isinstance(features[0], guppy.Multipoint):
            writer.shapeType = shapefile.MULTIPOINT
        else:
            raise NotImplementedError("cannot save type "
                                      "{0}".format(type(feature)))

        for feature in features:
            if writer.shapeType == shapefile.POLYLINE:
                parts = list_parts(feature)
            else:
                parts = [feature.vertices]
            writer.poly(parts)

    try:
        files = open_file_dict(fnms, 'wb')
        writer.saveShp(files['shp'])
        writer.saveShx(files['shx'])
        writer.saveDbf(files['dbf'])

    except Exception as e:
        raise e

    finally:
        for f in files.values():
            f.close()

