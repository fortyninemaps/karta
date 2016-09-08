"""
Data input functions for Karta

Used available drivers to read input data and return Karta geometry objects.
"""
import os
from numbers import Number
from osgeo import ogr
from . import geometry
from . import geojson
from . import shp
from . import gpx
from . import xyfile
from ..crs import GeographicalCRS, ProjectedCRS, LonLatWGS84
from .. import errors

### GeoInterface functions ###

def from_shape(obj, properties=None):
    """ Read a __geo_interface__ dictionary and return an appropriate karta
    object """
    return _from_shape(obj.__geo_interface__, None)

def _from_shape(d, properties):
    if d is None:
        return None
    elif d["type"] == "Feature":
        p = d["properties"]
        return _from_shape(d["geometry"], p)
    else:
        c = d["coordinates"]
        if c is None:
            return None
        elif d["type"] == "Point":
            return geometry.Point(c, properties=properties)
        elif d["type"] == "LineString":
            return geometry.Line(c, properties=properties)
        elif d["type"] == "Polygon":
            subs = [geometry.Polygon(c[i]) for i in range(1, len(c))]
            return geometry.Polygon(c[0], properties=properties, subs=subs)
        elif d["type"] == "MultiPoint":
            return geometry.Multipoint(c, properties=properties)
        elif d["type"] == "MultiLineString":
            return geometry.Multiline(c, properties=properties)
        elif d["type"] == "MultiPolygon":
            return geometry.Multipolygon(c, properties=properties)
        else:
            raise NotImplementedError("Geometry type {0} not "
                                      "implemented".format(d["type"]))

### GeoJSON functions ###

def read_geojson(f, crs=LonLatWGS84):
    """ Parse GeoJSON and return a list of geometries.

    f : file-like object or str
        file object to read from or a GeoJSON string
    crs : karta.crs.CRS
        CRS object to bind to new geometries
    """
    def _convert_crs(crsdict):
        # Deprecated
        if crsdict.get("type", None) not in ("name", "link"):
            crs = LonLatWGS84
        elif crsdict["type"] == "name":
            crs = geojson.GeoJSONNamedCRS(crsdict["properties"]["name"])
        elif crsdict["type"] == "link":
            crs = geojson.GeoJSONLinkedCRS(crsdict["properties"]["href"],
                                           crsdict["properties"]["type"])
        return crs

    def convert(geom, **kw):
        if isinstance(geom, geojson.Feature):
            res = [convert_feature(geom, **kw)]
        elif isinstance(geom, geojson.FeatureCollection):
            res = [convert_feature(f, **kw) for f in geom.features]
        elif isinstance(geom, geojson.GeometryCollection):
            res = [convert(item, **kw) for item in geom.geometries]
        else:
            res = convert_geometry(geom, **kw)
        return res

    def convert_geometry(geom, **kw):
        # TODO: Requires clean-up; was written before properties vs data was
        # clearly worked out. Hacks to bring it up to speed are ugly.
        if isinstance(geom, geojson.Point):
            kw.setdefault("properties", {})
            kw.setdefault("data", {})
            kw["properties"].update(kw["data"])
            del kw["data"]
            return geometry.Point(geom.coordinates, **kw)
        elif isinstance(geom, geojson.LineString):
            kw.setdefault("properties", {})
            kw.setdefault("data", {})
            kw["properties"].update(kw["data"])
            del kw["data"]
            return geometry.Line(geom.coordinates, **kw)
        elif isinstance(geom, geojson.Polygon):
            kw.setdefault("properties", {})
            kw.setdefault("data", {})
            kw["properties"].update(kw["data"])
            del kw["data"]
            return geometry.Polygon(geom.coordinates[0],
                                    subs=geom.coordinates[1:],
                                    **kw)
        elif isinstance(geom, geojson.MultiPoint):
            return geometry.Multipoint(geom.coordinates, **kw)
        elif isinstance(geom, geojson.MultiLineString):
            return [geometry.Line(coords, **kw)
                    for coords in geom.coordinates]
        elif isinstance(geom, geojson.MultiPolygon):
            return [geometry.Polygon(coords[0], subs=coords[1:], **kw)
                    for coords in geom.coordinates]
        else:
            raise TypeError("{0} is a not a GeoJSON entity".format(type(geom)))

    def convert_feature(feat, **kw):
        data = feat.properties["vector"]
        if len(data) == 0:
            data = {}
        else:
            for key, val in data.items():
                if any(isinstance(a, Number) or hasattr(a, "dtype") for a in val):
                    for i in range(len(val)):
                        if val[i] is None:
                            val[i] = float('nan')
        prop = feat.properties["scalar"]
        kw["data"] = data
        kw["properties"] = prop
        return convert(feat.geometry, **kw)

    R = geojson.GeoJSONReader(f)
    geom = R.parse()
    return convert(geom, crs=crs)

def _geojson_properties2karta(properties, n):
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

def read_xyfile(f, delimiter='', header_rows=0, astype=geometry.Multipoint, coordrank=2):
    """ Read an ASCII delimited table and return a geometry object given by *astype*.
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
    if stem.endswith(".shp"):
        stem = stem[:-4]
    shp = stem + '.shp'
    shx = stem + '.shx'
    dbf = stem + '.dbf'
    if check:
        for fnm in (shp, shx, dbf):
            if not os.path.isfile(fnm):
                raise Exception('missing {0}'.format(fnm))
    return {'shp':shp, 'shx':shx, 'dbf':dbf}

def ogr_read_shapefile(stem):
    fnms = get_filenames(stem)
    driver = ogr.GetDriverByName("ESRI Shapefile")
    ds = driver.Open(fnms["shp"], 0)
    layer = ds.GetLayer()

    try:
        _geoms = [_from_shape(gi, p)
                 for (gi,p) in zip(shp.ogr_read_geometries(layer),
                                   shp.ogr_read_attributes(layer))]
        crs = ogr_parse_srs(layer)
    finally:
        del ds, driver, layer

    geoms = []
    for g in _geoms:
        if isinstance(g, list):
            for part in g:
                part.crs = crs
            geoms.append(g)
        elif g is not None:
            g.crs = crs
            geoms.append(g)
    return geoms

def ogr_parse_srs(lyr):
    """ Given an OGR type with a `GetSpatialRef` method, return a matching CRS
    object. """
    srs = lyr.GetSpatialRef()
    if srs is None:
        crs = LonLatWGS84
    else:
        name = srs.GetAttrValue('PROJCS')
        if srs.IsGeographic():
            spheroid = "+a={a} +f={f}".format(a=srs.GetSemiMajor(),
                                              f=1.0/srs.GetInvFlattening())
            crs = GeographicalCRS(spheroid, name)
        else:
            crs = ProjectedCRS(srs.ExportToProj4(), name=name)
    return crs

# convenience binding
read_shapefile = ogr_read_shapefile


### GPX functions ###

def read_gpx_waypts(fnm):
    gpx_doc = gpx.GPX(fnm)
    return [_waypt2pt(pt) for pt in gpx_doc.waypts]

def read_gpx_tracks(fnm):
    gpx_doc = gpx.GPX(fnm)
    return [_track2lines(trk) for trk in gpx_doc.tracks]

def _waypt2pt(waypt):
    return geometry.Point(waypt.lonlat,
                          properties=waypt.properties,
                          crs=LonLatWGS84)

def _seg2line(seg):
    return geometry.Line([pt.lonlat for pt in seg.trkpts],
                         properties=seg.properties,
                         crs=LonLatWGS84)

def _track2lines(track):
    return [_seg2line(seg) for seg in track.trksegs]

