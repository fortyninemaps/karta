""" Provides functions for reading and writing ESRI shapefiles and returning a
geometry object. """

import os
import datetime
import numbers
from .. import errors

try:
    import osgeo
    from osgeo import ogr
    HAS_OSGEO = True
except ImportError:
    HAS_OSGEO = False

rankerror = IOError("feature must be in two or three dimensions to write "
                    "as shapefile")

if HAS_OSGEO:
    GEOMTYPES = {}
    for k in ogr.__dict__:
        if k.startswith("wkb"):
            GEOMTYPES[ogr.__dict__[k]] = k[3:]

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

def ogr_get_fieldtype(val):
    """ Return an ogr type integer that describes the type of *a*. """
    if isinstance(val, numbers.Number) or _isnumpytype(val):
        if isinstance(val, numbers.Integral) or _isnumpyint(val):
            return 0
        elif isinstance(val, numbers.Real) or _isnumpyfloat(val):
            return 2
        else:
            raise TypeError("cannot choose the correct dBase type for "
                            "{0}\n".format(type(value)))
    elif isinstance(val, str):
        return 6
    elif isinstance(val, bytes):
        return 8
    elif isinstance(val, datetime.datetime):
        return 11
    elif isinstance(val, datetime.date):
        return 9
    elif isinstance(val, datetime.time):
        return 10
    elif isinstance(val, list):
        try:
            t = ogr_get_fieldtype(val[0])
            if t == 0:
                return 1
            elif t == 2:
                return 3
            elif t == 6:
                return 7
            else:
                raise TypeError
        except TypeError as e:
            raise TypeError("cannot choose the correct dBase type for a list "
                            "with elements of type{0}".format(type(val[0])))
    else:
        raise TypeError("cannot choose the correct dBase type for "
                        "{0}".format(type(val)))

def get_geometry_type(gi):
    """ Return the geometry type from a __geo_interface__ dictionary """
    if gi["type"] == "Feature":
        return get_geometry_type(gi["geometry"])
    elif gi["type"] in ("FeatureCollection", "GeometryCollection"):
        return get_geometry_type(gi["geometries"][0])
    else:
        return gi["type"]

def ogr_get_polygon_points(poly):
    geoms = [poly.GetGeometryRef(i) for i in range(poly.GetGeometryCount())]
    pts = [geom.GetPoints()[:-1] for geom in geoms]
    return pts

def ogr_read_geometry(geom):
    # convert this to ogr_read_singlepart and write and additional ogr_read_multipart
    wkbtype = GEOMTYPES[geom.GetGeometryType()]
    if wkbtype in ('Point', 'Point25D', 'NDR', 'XDR'):
        jsontype = 'Point'
        pts = geom.GetPoint()
    elif wkbtype in ('LineString', 'LineString25D'):
        jsontype = 'LineString'
        pts = geom.GetPoints()
    elif wkbtype in ('Polygon', 'Polygon25D'):
        jsontype = 'Polygon'
        pts = ogr_get_polygon_points(geom)
    elif wkbtype in ('MultiPoint', 'MultiPoint25D'):
        jsontype = 'MultiPoint'
        pts = geom.GetPoints()
    elif wkbtype in ('MultiLineString', 'MultiLineString25D'):
        jsontype = 'MultiLineString'
        pts = [geom.GetGeometryRef(i).GetPoints()
               for i in range(geom.GetGeometryCount())]
    elif wkbtype in ('MultiPolygon', 'MultiPolygon25D'):
        jsontype = 'MultiPolygon'
        pts = [ogr_get_polygon_points(geom.GetGeometryRef(i))
               for i in range(geom.GetGeometryCount())]
    else:
        raise TypeError("WKB type {0} not handled".format(wkbtype))

    xmin, xmax, ymin, ymax = geom.GetEnvelope()
    return {'type': jsontype,
            'bbox': (xmin, ymin, xmax, ymax),
            'coordinates': pts}

def ogr_read_geometries(lyr):
    """ Given an ogr.Layer instance, return a list of geointerface dictionaries
    """
    gi_dicts = []
    for feature in lyr:
        geom = feature.GetGeometryRef()
        gi_dicts.append(ogr_read_geometry(geom))
    return gi_dicts

def ogr_read_attributes(lyr):
    """ Given an ogr.Layer instance, return a list of dictionaries with feature
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
            try:
                p[fn] = t(lyr.GetFeature(i).GetField(j))
            except TypeError:
                p[fn] = None
            except UnicodeDecodeError:
                # This is due to GDAL bug #5811 on Python 3
                # https://trac.osgeo.org/gdal/ticket/5811
                p[fn] = "gdal bug 5811"
        propertylist.append(p)

    return propertylist

def ogr_write(fnm, *objs):
    """ Write features to shapefile using OGR backend. Features may be karta
    geometry objects or __geo_interface__ features. """
    objs = [obj.__geo_interface__ if hasattr(obj, "_geotype") else obj
            for obj in objs]

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
    if len(proj4) == 0:
        proj4 = "+proj=longlat +datum=WGS84"
    srs = osgeo.osr.SpatialReference()
    try:
        srs.ImportFromProj4(proj4)
    except RuntimeError:
        raise errors.CRSError("invalid projection string: '{0}'".format(proj4))
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

write_shapefile = ogr_write
