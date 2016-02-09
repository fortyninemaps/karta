""" Provides functions for reading and writing ESRI shapefiles and returning a
geometry object. """

import os
import sys
import datetime
import numbers
import shapefile
from .. import errors

try:
    import osgeo
    from osgeo import ogr
    HAS_OSGEO = True
except ImportError:
    HAS_OSGEO = False

rankerror = IOError("feature must be in two or three dimensions to write "
                    "as shapefile")

GEOMTYPES = {1: "Point",
             2: "LineString",
             3: "Polygon",
             4: "MultiPoint",
             5: "MultiLineString",
             6: "MultiPolygon"}

OGRTYPES = {"Point": 1,
            "LineString": 2,
            "Polygon": 3,
            "MultiPoint": 4,
            "MultiLineString": 5,
            "MultiPolygon": 6}

FIELDTYPES = {0: int,
              1: list,      # IntegerList
              2: float,
              3: list,      # RealList
              4: str,
              5: list,      # StringList
              6: str,       # WideString
              7: list,      # WideStringList
              8: bin,
              9: datetime.date,
              10: datetime.time,
              11: datetime.datetime}

OGRFIELDTYPES = {int: 0,
                 float: 2,
                 str: 4,
                 bin: 8,
                 datetime.date: 9,
                 datetime.time: 10,
                 datetime.datetime: 11}

def ogr_get_fieldtype(a):
    """ Return an ogr type integer that describes the type of *a*. """
    if isinstance(a, numbers.Number) or _isnumpytype(a):
        if isinstance(a, numbers.Integral) or _isnumpyint(a):
            return 0
        elif isinstance(a, numbers.Real) or _isnumpyfloat(a):
            return 2
        else:
            raise TypeError("cannot choose the correct dBase type for "
                            "{0}\n".format(type(value)))
    elif isinstance(a, str):
        return 4
    elif isinstance(a, datetime.date):
        return 9
    elif isinstance(a, datetime.time):
        return 10
    elif isinstance(a, datetime.datetime):
        return 11
    else:
        raise TypeError("cannot choose the correct dBase type for "
                        "{0}\n".format(type(value)))

def get_geometry_type(gi):
    """ Return the geometry type from a __geo_interface__ dictionary """
    if gi["type"] == "Feature":
        return get_geometry_type(gi["geometry"])
    elif gi["type"] in ("FeatureCollection", "GeometryCollection"):
        return get_geometry_type(gi["geometries"][0])
    else:
        return gi["type"]

def ogr_read_geometries(lyr):
    """ Given an ogr.Layer instance, return a list of geointerace dictionaries """
    gi_dicts = []

    for feature in lyr:
        geom = feature.GetGeometryRef()
        t = GEOMTYPES[geom.GetGeometryType()]
        if t == 'Point':
            pts = geom.GetPoint()
        elif t == 'Polygon':
            nrings = geom.GetGeometryCount()
            outer_ring = geom.GetGeometryRef(0)
            if nrings != 1:
                inner_rings = [geom.GetGeometryRef(i) for i in range(1, nrings)]
            else:
                inner_rings = []
            pts = [outer_ring.GetPoints()[:-1]]
            for ring in inner_rings:
                pts.append(ring.GetPoints()[:-1])
        else:
            pts = geom.GetPoints()

        xmin, xmax, ymin, ymax = geom.GetEnvelope()
        d = {'type': t,
             'bbox': (xmin, ymin, xmax, ymax),
             'coordinates': pts}
        gi_dicts.append(d)
    return gi_dicts

def ogr_read_attributes(lyr):
    """ Given an ogr.Layer instance, return a list of dictionaries with fetaure
    property data. """
    lyr_defn = lyr.GetLayerDefn()
    nfields = lyr_defn.GetFieldCount()
    field_defns = [lyr_defn.GetFieldDefn(i) for i in range(nfields)]
    fieldnames = [fd.GetName() for fd in field_defns]
    fieldtypes = [FIELDTYPES[fd.GetType()] for fd in field_defns]

    nfeat = lyr.GetFeatureCount()
    propertylist = []
    for i in range(nfeat):
        p = {}
        for j, (fn, t) in enumerate(zip(fieldnames, fieldtypes)):
            p[fn] = t(lyr.GetFeature(i).GetField(j))
        propertylist.append(p)

    return propertylist

def ogr_read_attribute_table(lyr):
    """ Given an ogr.Layer instance, return a list of field names, a list of
    field types, and a list of field data lists. """
    lyr_defn = lyr.GetLayerDefn()
    nfields = lyr_defn.GetFieldCount()
    field_defns = [lyr_defn.GetFieldDefn(i) for i in range(nfields)]
    fieldnames = [fd.GetName() for fd in field_defns]
    fieldtypes = [FIELDTYPES[fd.GetType()] for fd in field_defns]

    nfeat = lyr.GetFeatureCount()
    fielddata = []
    for i in range(nfeat):
        feat = lyr.GetFeature(i)
        fielddata.append([t(feat.GetField(j)) for (t, j)
                          in zip(fieldtypes, range(nfields))])

    return fieldnames, fieldtypes, fielddata

def ogr_write(fnm, *objs):
    driver = ogr.GetDriverByName("ESRI Shapefile")
    if driver is None:
        raise ValueError("failure loading OGR driver for 'ESRI Shapefile'")
    if os.path.isfile(fnm):
        driver.DeleteDataSource(fnm)

    ds = driver.CreateDataSource(fnm)
    if ds is None:
        raise ValueError("failure creating ogr.DataSource from ogr.Driver")

    # create a spatial reference
    proj4 = objs[0]["properties"]["_karta_proj4"]\
                .replace("lonlat", "longlat")\
                .replace("latlon", "latlong")
    srs = osgeo.osr.SpatialReference()
    srs.ImportFromProj4(proj4)
    if srs is None:
        raise ValueError("failure creating osr.SpatialReference")

    # create layer
    geo_t = get_geometry_type(objs[0])
    ogr_t = OGRTYPES[geo_t]
    layer_name = os.path.splitext(os.path.split(fnm)[1])[0]
    lyr = ds.CreateLayer(layer_name, srs, ogr_t)
    if lyr is None:
        raise ValueError("failure creating ogr.Layer within ogr.DataSource")

    # add attribute table, reserving the first column for 'id'
    idfield = ogr.FieldDefn("id", ogr.OFTInteger)
    lyr.CreateField(idfield)

    keys = set.intersection(*[set(obj["properties"].keys()) for obj in objs])
    keys = [k for k in keys if not k.startswith("_karta")]
    for k in keys:
        t = ogr_get_fieldtype(obj["properties"][k][0])
        lyr.CreateField(ogr.FieldDefn(k, t))

    # add geometries
    for i, obj in enumerate(objs):
        ogr_write_feature(lyr, obj, id=i)
    return

def ogr_write_feature(lyr, gi, id=0):
    """ Write the geometry encoded in a __geointerface__ dictionary to a
    shapefile. """
    feature = ogr.Feature(lyr.GetLayerDefn())
    feature.SetField("id", id)
    for i in range(1, lyr.GetLayerDefn().GetFieldCount()):
        attr_name = lyr.GetFieldDefn(i).GetNameRef()
        feature.SetField(attr_name, gi["properties"][attr_name])

    if gi["geometry"]["type"] == "Polygon":
        ogr_write_ring_geometry(feature, gi["geometry"])
    else:
        ogr_write_geometry(feature, gi["geometry"])
    lyr.CreateFeature(feature)
    return

def ogr_write_geometry(feature, gi):
    geom = ogr.Geometry(OGRTYPES[gi["type"]])
    for pt in gi["coordinates"]:
        geom.AddPoint(pt[0], pt[1])
    feature.SetGeometry(geom)
    return

def ogr_write_ring_geometry(feature, gi):
    geom = ogr.Geometry(OGRTYPES[gi["type"]])
    for poly in gi["coordinates"]:
        ring = ogr.Geometry(ogr.wkbLinearRing)
        for pt in poly:
            ring.AddPoint(pt[0], pt[1])
        geom.AddGeometry(ring)
    geom.CloseRings()
    feature.SetGeometry(geom)
    return

def _isnumpyint(o):
    return hasattr(o, "dtype") and \
            o.dtype in ("int8", "int16", "int32", "int64")

def _isnumpyfloat(o):
    return hasattr(o, "dtype") and \
            o.dtype in ("float16", "float32", "float64", "float128")

def _isnumpytype(o):
    return hasattr(o, "dtype")

### connections for the old pyshp backend ###

import shapefile

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
    """ Add geometry properties *properties* to shapefile writer *writer*. """
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
    """ Add *points* data to shapefile writer *writer*. """
    writer.field("ID", "I", "8")
    if points.data is not None:
        keys = points.data.fields
        for key in keys:
            t = property_field_type(points.data.get(0)[key])
            dec = 0 if t not in ("N", "F", "O", "I") else 16
            writer.field(key, fieldType=t, size="100", decimal=dec)
        for i, pt in enumerate(points):
            writer.record(str(i), *[pt.data.get(0)[key] for key in keys])
    else:
        for i in range(len(points)):
            writer.record(str(i))
    return

def write_multipoint2(mp, fstem):
    w = shapefile.Writer(shapeType=shapefile.POINT)
    for vertex in mp.vertices:
        w.point(*vertex)
    addfields_points(w, mp)
    w.save(fstem)
    write_projection(mp, fstem)
    return

def write_multipoint3(mp, fstem):
    w = shapefile.Writer(shapeType=shapefile.POINTZ)
    for vertex in mp.vertices:
        w.point(*vertex)
    addfields_points(w, mp)
    w.save(fstem)
    write_projection(mp, fstem)
    return

def write_line2(line, fstem):
    w = shapefile.Writer(shapeType=shapefile.POLYLINE)
    w.poly(shapeType=shapefile.POLYLINE, parts=[line.vertices])
    addfields(w, line.properties)
    w.save(fstem)
    write_projection(line, fstem)
    return

def write_line3(line, fstem):
    w = shapefile.Writer(shapeType=shapefile.POLYLINEZ)
    w.poly(shapeType=shapefile.POLYLINEZ, parts=[line.vertices])
    addfields(w, line.properties)
    w.save(fstem)
    write_projection(line, fstem)
    return

def write_poly2(poly, fstem):
    w = shapefile.Writer(shapeType=shapefile.POLYGON)
    w.poly(shapeType=shapefile.POLYGON, parts=[[v for v in poly.vertices]])
    addfields(w, poly.properties)
    w.save(fstem)
    write_projection(poly, fstem)
    return

def write_poly3(poly, fstem):
    w = shapefile.Writer(shapeType=shapefile.POLYGONZ)
    w.poly(shapeType=shapefile.POLYGONZ, parts=[[v for v in poly.vertices]])
    addfields(w, poly.properties)
    w.save(fstem)
    write_projection(poly, fstem)
    return

def shapefile_type(feature):
    """ Match a feature with a shapefile type identifier """
    if feature._geotype == "Multipoint":
        if feature.rank == 2:
            return shapefile.POINT
        elif feature.rank == 3:
            return shapefile.POINTZ
        else:
            raise rankerror
    elif feature._geotype == "Line":
        if feature.rank == 2:
            return shapefile.POLYLINE
        elif feature.rank == 3:
            return shapefile.POLYLINEZ
        else:
            raise rankerror
    elif feature._geotype == "Polygon":
        if feature.rank == 2:
            return shapefile.POLYGON
        elif feature.rank == 3:
            return shapefile.POLYGONZ
        else:
            raise rankerror
    else:
        raise TypeError("Expected list with items of geotype 'Multipoint',"
                        "'Line', or 'Polygon'")

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

    shapeType = shapefile_type(features[0])
    w = shapefile.Writer(shapeType=shapeType)
    if shapeType in (shapefile.POINT, shapefile.POINTZ):

        # add geometry
        for feature in features:
            for pt in feature:
                w.point(*pt.vertex)

        # add records
        w.field("ID", "I", "8")
        if features[0].data is not None:
            keys = set(features[0].data.fields)     # for testing similarity
            keylist = features[0].data.fields       # preserves order
        else:
            keys = []

        if len(keys) != 0 and all(keys == set(f.data.fields) for f in features[1:]):
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
            w.poly([feature.vertices])

        # add records
        w.field("ID", "I", "8")
        keys = set(features[0].properties.keys())   # for testing similarity
        keylist = list(keys)                        # preserves order

        if len(keys) != 0 and all(keys == set(f.properties.keys()) for f in features[1:]):
            for key in features[0].properties:
                value = features[0].properties[key]
                propertyType = property_field_type(value)
                length = "100"
                w.field(key.upper(), propertyType, length)
            for (i, feature) in enumerate(features):
                values = []
                for key in keylist:
                    values.append(feature.properties[key])
                w.record(str(i), *values)
        else:
            for (i, feature) in enumerate(features):
                w.record(str(i))

    w.save(fstem)
    write_projection(features[0], fstem)
    return

def write_projection(geom, fstem):
    try:
        wkt = geom.crs.get_wkt()
        with open(fstem+".prj", "w") as f:
            f.write(wkt)
    except errors.CRSError:
        sys.stderr.write("Projection information not saved (requires osgeo)\n")
    except AttributeError as e:
        sys.stderr.write(str(e))
        sys.stderr.write("\n")

