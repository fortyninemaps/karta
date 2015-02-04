""" Provides functions for reading and writing ESRI shapefiles and returning a
geometry object. """

import datetime
import numbers
import shapefile

def _isnumpyint(o):
    return hasattr(o, "dtype") and \
            o.dtype in ("int8", "int16", "int32", "int64")

def _isnumpyfloat(o):
    return hasattr(o, "dtype") and \
            o.dtype in ("float16", "float32", "float64", "float128")

def _isnumpytype(o):
    return hasattr(o, "dtype")

def property_field_type(value):
    """ Determine the appropriate dBase field type for *value* """
    if isinstance(value, numbers.Number) or _isnumpytype(value):
        ## pyshp doesn't handle the full variety of numeric types
        ## (float, double, long), only 'N'
        #if isinstance(value, numbers.Integral) or _isnumpyint(value):
        #    desc = "I"
        #elif isinstance(value, numbers.Real) or _isnumpyfloat(value):
        #    desc = "O"
        if isinstance(value, (numbers.Integral, numbers.Real)) \
                or _isnumpyint(value) or _isnumpyfloat(value):
            desc = "N"
        else:
            raise TypeError("cannot choose the correct dBase type for "
                            "{0}\n".format(type(value)))
    elif isinstance(value, str):
        desc = "C"
    elif isinstance(value, datetime.datetime):
        desc = "@"
    elif isinstance(value, datetime.date):
        desc = "D"
    elif isinstance(value, bool):
        desc = "L"
    else:
        raise TypeError("cannot choose the correct dBase type for "
                        "{0}\n".format(type(value)))
    return desc

def addfields(writer, properties):
    writer.field("ID", "I", "8")
    values = []
    for key in properties:
        value = properties[key]
        typ = property_field_type(value)
        dec = 0 if typ not in ("N", "F", "O", "I") else 16
        writer.field(key.upper(), fieldType=typ, size="100", decimal=dec)
        values.append(value)
    writer.record("0", *values)
    return

def addfields_points(writer, points):
    writer.field("ID", "I", "8")
    keys = list(points.data.keys())
    for key in keys:
        t = property_field_type(points.data[key][0])
        dec = 0 if t not in ("N", "F", "O", "I") else 16
        writer.field(key, fieldType=t, size="100", decimal=dec)
    for i, pt in enumerate(points):
        writer.record(str(i), *[pt.data[key] for key in keys])
    return

def write_multipoint2(mp, fstem):
    w = shapefile.Writer(shapeType=shapefile.POINT)
    for vertex in mp.vertices:
        w.point(*vertex)
    addfields_points(w, mp)
    w.save(fstem)
    return

def write_multipoint3(mp, fstem):
    w = shapefile.Writer(shapeType=shapefile.POINTZ)
    for vertex in mp.vertices:
        w.point(*vertex)
    addfields_points(w, mp)
    w.save(fstem)
    return

def write_line2(line, fstem):
    w = shapefile.Writer(shapeType=shapefile.POLYLINE)
    w.poly(shapeType=shapefile.POLYLINE, parts=[line.vertices])
    addfields(w, line.properties)
    w.save(fstem)
    return

def write_line3(line, fstem):
    w = shapefile.Writer(shapeType=shapefile.POLYLINEZ)
    w.poly(shapeType=shapefile.POLYLINEZ, parts=[line.vertices])
    addfields(w, line.properties)
    w.save(fstem)
    return

def write_poly2(poly, fstem):
    w = shapefile.Writer(shapeType=shapefile.POLYGON)
    w.poly(shapeType=shapefile.POLYGON, parts=[poly.vertices])
    addfields(w, poly.properties)
    w.save(fstem)
    return

def write_poly3(poly, fstem):
    w = shapefile.Writer(shapeType=shapefile.POLYGONZ)
    w.poly(shapeType=shapefile.POLYGONZ, parts=[poly.vertices])
    addfields(w, poly.properties)
    w.save(fstem)
    return

def write_shapefile(features, fstem):
    """ Write *features* to a shapefile. All features must be of the same type
    (i.e. Multipoints, Lines, Polygons). """
    if len(features) == 1:
        features[0].to_shapefile(fstem)
        return
    elif not all(f._geotype == features[0]._geotype for f in features[1:]):
        raise IOError("all features must be of the same type")
    elif not all(f.rank == features[0].rank for f in features[1:]):
        raise IOError("all features must have the same dimensionality")

    RankError = IOError("feature must be in two or three dimensions to write "
                        "as shapefile")
    if features[0]._geotype == "Multipoint":
        if features[0].rank == 2:
            typ = shapefile.POINT
        elif features[0].rank == 3:
            typ = shapefile.POINTZ
        else:
            raise RankError
    elif features[0]._geotype == "Line":
        if features[0].rank == 2:
            typ = shapefile.POLYLINE
        elif features[0].rank == 3:
            typ = shapefile.POLYLINEZ
        else:
            raise RankError
    elif features[0]._geotype == "Polygon":
        if features[0].rank == 2:
            typ = shapefile.POLYGON
        elif features[0].rank == 3:
            typ = shapefile.POLYGONZ
        else:
            raise RankError

    w = shapefile.Writer(shapeType=typ)
    if typ in (shapefile.POINT, shapefile.POINTZ):

        # add geometry
        for feature in features:
            for pt in feature:
                w.point(*pt.vertex)

        # add records
        w.field("ID", "I", "8")
        keys = set(features[0].data.keys())     # for testing similarity
        keylist = list(keys)                    # preserves order

        if len(keys) != 0 and all(keys == set(f.data.keys()) for f in features[1:]):
            for key in keylist:
                testvalue = features[0].data[key][0]
                w.field(key.upper(), property_field_type(testvalue), "100")
            i = 0
            for feature in features:
                for pt in feature:
                    w.record(str(i), *[pt.data[key] for key in keylist])
                    i += 1
        else:
            i = 0
            for feature in features:
                for pt in feature:
                    w.record(str(i))
                    i += 1

    else:

        # add geometry
        for feature in features:
            #print(feature)
            w.poly([feature.vertices])

        # add records
        w.field("ID", "I", "8")
        keys = set(features[0].properties.keys())   # for testing similarity
        keylist = list(keys)                        # preserves order

        if len(keys) != 0 and all(keys == set(f.data.keys()) for f in features[1:]):
            for key in features[0].properties:
                value = features[0].properties[key]
                typ = property_field_type(value)
                length = "100"
                w.field(key.upper(), typ, length)
            for (i, feature) in enumerate(features):
                values = []
                for key in keylist:
                    value.append(feature.properties[key])
                w.record(str(i), *values)
        else:
            for (i, feature) in enumerate(features):
                w.record(str(i))

    w.save(fstem)
    return

