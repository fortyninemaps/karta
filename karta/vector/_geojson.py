import picogeojson
from .utilities import _as_nested_lists, _reproject_nested
from ..crs import LonLatWGS84

def crs_from_urn(urn):
    return {'type': 'name', 'properties': {'name': urn}}

def crs_from_karta(crs):
    if hasattr(crs, "jsonhref") and hasattr(crs, "jsontype"):
        return {'type': 'link', 'properties': {'href': crs.jsonname, 'type': crs.jsontype}}
    elif hasattr(crs, "jsonname"):
        return {'type': 'name', 'properties': {'name': crs.jsonname}}
    else:
        # NOTE: make sure CRS gets written to file by caller
        return {'type': 'link', 'properties': {'href': "", 'type': 'proj4'}}

class GeoJSONOutMixin(object):
    """ Mixin class to be added to geometry objects, adding geojson
    functionality.
    """

    _serializer = picogeojson.Serializer(antimeridian_cutting=False,
                                         enforce_poly_winding=False)

    @staticmethod
    def _as_named_tuple(geom, **kwargs):
        """ Convert one or more Geometry instances to GeoJSON-structured named tuples.
        Parameters
        ----------
        *geoms : subtypes of karta.Geometry
            karta vector geometries to convert
        urn : str, optional
            URN defining a specific CRS to use
        Returns
        -------
        Either a Feature or a FeatureCollection
        Raises
        ------
        TypeError
            if one or more of geoms has an unrecognized `_geotype` attribute
        """
        if kwargs.get("urn", None) is not None:
            crs = crs_from_urn(kwargs["urn"])
        else:
            crs = crs_from_karta(geoms[0].crs)

        if geom._geotype == "Point":
            g = picogeojson.Point(geom.vertex, crs)
        elif geom._geotype == "Line":
            g = picogeojson.LineString(_as_nested_lists(geom.vertices), crs)
        elif geom._geotype == "Polygon":
            verts = [_as_nested_lists(geom.vertices_ring)]
            for sub in geom.subs:
                verts.append(_as_nested_lists(sub.vertices))
            g = picogeojson.Polygon(verts, crs)
        elif geom._geotype == "Multipoint":
            g = picogeojson.MultiPoint(_as_nested_lists(geom.vertices), crs)
        elif geom._geotype == "Multiline":
            g = picogeojson.MultiLineString(_as_nested_lists(geom.vertices), crs)
        elif geom._geotype == "Multipolygon":
            g = picogeojson.MultiPolygon(_as_nested_lists(geom.vertices_ring), crs)
        else:
            raise TypeError("unhandled type: {0}".format(type(geom)))

        data = {}
        if hasattr(geom, "data"):
            for field in geom.data.fields:
                data[field] = geom.d[field]

        properties = geom.properties
        properties.update(data)
        return picogeojson.Feature(g, properties)

    def as_geojson(self, indent=None, urn=None, force_wgs84=True):
        """ Output representation of internal data as a GeoJSON string.
        Parameters
        ----------
        indent : int, optional
            indentation of generated GeoJSON (default 2)
        force_wgs84 : bool, optional
            Forces output to use geographical coordinates with the WGS84 datum,
            as recommended by the GeoJSON draft specification
            (https://datatracker.ietf.org/doc/draft-ietf-geojson/).
            If *urn* is not set, "urn:ogc:def:crs:OGC:1.3:CRS84" is used.
            (default True)
        urn : str, optional
            overrides GeoJSON CRS with provided URN string
        """
        if force_wgs84 and (urn is None):
            urn = "urn:ogc:def:crs:OGC:1.3:CRS84"
        if force_wgs84 and (self.crs != LonLatWGS84):
            kw = dict(properties=self.properties, crs=LonLatWGS84)
            if hasattr(self, "data"):
                kw["data"] = self.data
            if hasattr(self, "vertices"):
                vertices = _reproject_nested(self.vertices, self.crs, LonLatWGS84)
            else:
                vertices = _reproject(self.vertex, self.crs, LonLatWGS84)
            geo = type(self)(vertices, **kw)
        else:
            geo = self

        return self._serializer(self._as_named_tuple(geo, urn=urn),
                                indent=indent)

    def to_geojson(self, f, indent=None, urn=None, force_wgs84=True):
        """ Write data as a GeoJSON string to a file-like object `f`.
        Parameters
        ----------
        f : str or file-like object
            file to receive GeoJSON string
        indent : int, optional
            indentation of generated GeoJSON (default None)
        force_wgs84 : bool, optional
            Forces output to use geographical coordinates with the WGS84 datum,
            as recommended by the GeoJSON draft specification
            (https://datatracker.ietf.org/doc/draft-ietf-geojson/).
            If *urn* is not set, "urn:ogc:def:crs:OGC:1.3:CRS84" is used.
            (default True)
        urn : str, optional
            overrides GeoJSON CRS with provided URN string
        """
        try:
            if not hasattr(f, "write"):
                fobj = open(f, "w")
            else:
                fobj = f
            fobj.write(self.as_geojson(indent=indent, urn=urn,
                                       force_wgs84=force_wgs84))
        finally:
            if not hasattr(f, "write"):
                fobj.close()
        return

