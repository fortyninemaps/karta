"""
Shapefile driver for Karta using OGR backend

Overview
--------

`ogr_write` takes geointerface dictionaries or Karta objects and writes to ESRI
shapefile

`ogr_read_geometry` converts an OGR geometry object to a geointerface dictionary
"""

import os
import datetime
import numbers
from .. import errors
import osgeo
from osgeo import ogr

OGRGEOMTYPES = {}
for k in ogr.__dict__:
    if k.startswith("wkb"):
        OGRGEOMTYPES[ogr.__dict__[k]] = k[3:]

OGRTYPES = {"Point": ogr.wkbPoint,
            "LineString": ogr.wkbLineString,
            "Polygon": ogr.wkbPolygon,
            "MultiPoint": ogr.wkbMultiPoint,
            "MultiLineString": ogr.wkbMultiLineString,
            "MultiPolygon": ogr.wkbMultiPolygon}

# dictionary that matches each OGR type to a Python builtin
FIELDTYPES = {ogr.OFTInteger: int,
              ogr.OFTIntegerList: list,
              ogr.OFTReal: float,
              ogr.OFTRealList: list,
              ogr.OFTString: str,
              ogr.OFTStringList: list,
              ogr.OFTWideString: str,
              ogr.OFTWideStringList: list,
              ogr.OFTBinary: bin,
              ogr.OFTDate: datetime.date,
              ogr.OFTTime: datetime.time,
              ogr.OFTDateTime: datetime.datetime}

try:
    FIELDTYPES[ogr.OFTInteger64] = int
    FIELDTYPES[ogr.OFTInteger64List] = list
except AttributeError:
    pass

def ogr_get_fieldtype(val):
    """ Returns a tuple of
        - an OGR type integer that describes the type of *a*
        - a default width that should work for most cases
    """
    if isinstance(val, numbers.Number) or _isnumpytype(val):
        if isinstance(val, numbers.Integral) or _isnumpyint(val):
            return ogr.OFTInteger, 32
        elif isinstance(val, numbers.Real) or _isnumpyfloat(val):
            return ogr.OFTReal, 32
        else:
            raise TypeError("cannot choose the correct dBase type for "
                            "{0}\n".format(type(val)))
    elif isinstance(val, str):
        return ogr.OFTString, 180
    elif isinstance(val, bytes):
        return ogr.OFTBinary, 32
    elif isinstance(val, datetime.datetime):
        return ogr.OFTDateTime, 64
    elif isinstance(val, datetime.date):
        return ogr.OFTDate, 32
    elif isinstance(val, datetime.time):
        return ogr.OFTTime, 32
    elif isinstance(val, list):
        t = ogr_get_fieldtype(val[0])[0]
        if t == ogr.OFTInteger:
            return ogr.OFTIntegerList, 1000
        elif t == ogr.OFTReal:
            return ogr.OFTRealList, 1000
        elif t == ogr.OFTString:
            return ogr.OFTStringList, 1000
        else:
            raise TypeError("cannot choose the correct dBase type for a list "
                            "with elements of type '{0}'".format(type(val[0])))
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
    """ Return a list of lists of ogr.Polygon coordinate rings. """
    geoms = [poly.GetGeometryRef(i) for i in range(poly.GetGeometryCount())]
    pts = [geom.GetPoints()[:-1] for geom in geoms]
    return pts

def ogr_read_geometry(geom):
    """ Convert an ogr.Geometry to a __geo_interface__ dictionary. """
    if geom is None:
        return None
    wkbtype = OGRGEOMTYPES[geom.GetGeometryType()]
    if wkbtype in ('Point', 'Point25D', 'PointZ', 'PointM', 'PointZM',
                   'NDR', 'XDR'):
        jsontype = 'Point'
        pts = geom.GetPoint()
    elif wkbtype in ('LineString', 'LineString25D', 'LineStringZ',
                     'LineStringM', 'LineStringZM'):
        jsontype = 'LineString'
        pts = geom.GetPoints()
    elif wkbtype in ('Polygon', 'Polygon25D', 'PolygonZ', 'PolygonM',
                     'PolygonZM'):
        jsontype = 'Polygon'
        pts = ogr_get_polygon_points(geom)
    elif wkbtype in ('MultiPoint', 'MultiPoint25D', 'MultiPointZ',
                     'MultiPointM', 'MultiPointZM'):
        jsontype = 'MultiPoint'
        pts = [geom.GetGeometryRef(i).GetPoints()[0]
               for i in range(geom.GetGeometryCount())]
    elif wkbtype in ('MultiLineString', 'MultiLineString25D', 'MultiLinStringZ', 
                     'MultiLineStringM', 'MultiLineStringZM'):
        jsontype = 'MultiLineString'
        pts = [geom.GetGeometryRef(i).GetPoints()
               for i in range(geom.GetGeometryCount())]
    elif wkbtype in ('MultiPolygon', 'MultiPolygon25D', 'MultiPolygonZ',
                     'MultiPolygonM', 'MultiPolygonZM'):
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
    """ Given an ogr.Layer instance, return a list of __geo_interface__
    dictionaries.
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
    """ Write features to shapefile using OGR backend. Features are
    __geo_interface__ Feature mappings. """
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
        typ, default_width = ogr_get_fieldtype(objs[0]["properties"][k])
        fd = ogr.FieldDefn(k, typ)
        fd.SetWidth(default_width)
        lyr.CreateField(fd)

    # add geometries
    for i, obj in enumerate(objs):
        ogr_write_feature(lyr, obj, id=i)
    return

def ogr_write_feature(lyr, gi, id=0):
    """ Write the geometry encoded in a __geointerface__ dictionary to an OGR
    Layer. """
    layer_def = lyr.GetLayerDefn()
    feature = ogr.Feature(layer_def)
    feature.SetField("id", id)
    for i in range(1, layer_def.GetFieldCount()):
        attr_name = layer_def.GetFieldDefn(i).GetNameRef()
        feature.SetField(attr_name, gi["properties"][attr_name])

    if gi["geometry"]["type"] == "Point":
        geom = ogr_asPoint(gi["geometry"])
    elif gi["geometry"]["type"] == "LineString":
        geom = ogr_asLineString(gi["geometry"])
    elif gi["geometry"]["type"] == "Polygon":
        geom = ogr_asPolygon(gi["geometry"])
    elif gi["geometry"]["type"] == "MultiPoint":
        geom = ogr_asMultiPoint(gi["geometry"])
    elif gi["geometry"]["type"] == "MultiLineString":
        geom = ogr_asMultiLineString(gi["geometry"])
    elif gi["geometry"]["type"] == "MultiPolygon":
        geom = ogr_asMultiPolygon(gi["geometry"])
    else:
        raise ValueError("geometry type unhandled: {0}".format(gi["geometry"]["type"]))
    feature.SetGeometry(geom)
    lyr.CreateFeature(feature)
    return

def ogr_asPoint(gi):
    geom = ogr.Geometry(1)
    geom.AddPoint(gi["coordinates"][0], gi["coordinates"][1])
    return geom

def ogr_asLineString(gi):
    geom = ogr.Geometry(2)
    for pt in gi["coordinates"]:
        geom.AddPoint(float(pt[0]), float(pt[1]))
    return geom

def ogr_asPolygon(gi):
    geom = ogr.Geometry(3)
    for ring in gi["coordinates"]:
        ogr_ring = ogr.Geometry(ogr.wkbLinearRing)
        for pt in ring:
            ogr_ring.AddPoint(pt[0], pt[1])
        geom.AddGeometry(ogr_ring)
    geom.CloseRings()
    return geom

def ogr_asMultiPoint(gi):
    geom = ogr.Geometry(4)
    for pt in gi["coordinates"]:
        g = ogr.Geometry(1)
        g.AddPoint(float(pt[0]), float(pt[1]))
        geom.AddGeometry(g)
    return geom

def ogr_asMultiLineString(gi):
    geom = ogr.Geometry(5)
    for linestring in gi["coordinates"]:
        ogr_linstring = ogr.Geometry(2)
        for pt in linestring:
            ogr_linstring.AddPoint(float(pt[0]), float(pt[1]))
        geom.AddGeometry(ogr_linstring)
    return geom

def ogr_asMultiPolygon(gi):
    geom = ogr.Geometry(5)
    for polygon in gi["coordinates"]:
        ogr_polygon = ogr.Geometry(3)
        for ring in polygon:
            ogr_ring = ogr.Geometry(ogr.wkbLinearRing)
            for pt in ring:
                ogr_ring.AddPoint(pt[0], pt[1])
            ogr_polygon.AddGeometry(ogr_ring)
        ogr_polygon.CloseRings()
        geom.AddGeometry(ogr_polygon)
    return geom

def _isnumpyint(o):
    return hasattr(o, "dtype") and \
            o.dtype in ("int8", "int16", "int32", "int64")

def _isnumpyfloat(o):
    return hasattr(o, "dtype") and \
            o.dtype in ("float16", "float32", "float64", "float128")

def _isnumpytype(o):
    return hasattr(o, "dtype")

write_shapefile = ogr_write

def write_shapefile(fnm, *geometries):
    return ogr_write(fnm, *[g.__geo_interface__ for g in geometries])

class ShapefileOutMixin(object):
    """ Mixin class to be added to geometry objects, adding shapefile
    functionality.
    """

    def to_shapefile(self, fnm):
        """ Save line to a shapefile """
        if not fnm.endswith(".shp"):
            fnm = fnm + ".shp"
        ogr_write(fnm, self.__geo_interface__)
        return

