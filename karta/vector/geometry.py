"""
Geographical measurement and analysis

Provides Point, Line, and Polygon classes, and their Multipart equivalents,
with methods for simple measurements such as distance, area, and direction.
"""
from __future__ import division

import math
import itertools
import numbers
import numpy as np
from .decorators import cache_decorator
from .geojson import GeoJSONOutMixin
from .shp import ShapefileOutMixin
from . import xyfile
from .table import Table, Indexer
from .utilities import _reproject, _reproject_nested, _flatten, _as_nested_lists
from .coordstring import CoordString
from .rtree import RTree
from .quadtree import QuadTree
from . import vectorgeo as _cvectorgeo
from . import dateline as _cdateline
from . import intersection as _cintersection
from . import convexhull as _cconvexhull
from . import contains as _ccontains
from .. import geodesy
from ..crs import Cartesian, CartesianCRS, GeographicalCRS, LonLatWGS84
from ..crs import SphericalEarth
from ..errors import GeometryError, GGeoError, CRSError

class Geometry(object):
    """ Abstract base class for all geometry types """

    #__slots__ = ["_geotype", "properties", "crs", "_cache"]

    def __init__(self, properties=None, crs=Cartesian):
        self.crs = crs
        if isinstance(properties, dict):
            self.properties = properties
        elif properties is not None:
            raise TypeError("properties must be a dictionary")
        else:
            self.properties = {}
        self._cache = {}
        self._geotype = None
        return

class Point(Geometry, GeoJSONOutMixin, ShapefileOutMixin):
    """ Point object instantiated with:

    Parameters
    ----------
    coords : 2-tuple or 3-tuple
    properties : dict, optional
        geometry specific metadata (default None)
    crs : karta.crs.CRS, optional
        coordinate system for geometry (default Cartesian)
    """
    #__slots__ = ["vertex"]

    def __init__(self, coords, properties=None, **kwargs):
        if len(coords) not in (2, 3):
            raise TypeError("Point coordinates must be a sequence")
        super(Point, self).__init__(properties=properties, **kwargs)
        self.vertex = tuple(coords)
        self._geotype = "Point"
        return

    def __getitem__(self, idx):
        return self.vertex[idx]

    def __repr__(self):
        return 'Point({0}, {1})'.format(self.x, self.y)

    def __eq__(self, other):
        try:
            return (tuple(self.vertex) == tuple(other.vertex)) and \
                   (self.properties == other.properties) and \
                   (self.crs == other.crs)
        except AttributeError:
            return False

    def __neq__(self, other):
        return ~(self == other)

    def __hash__(self):
        ht = hash(self._geotype)
        hd = hash(self.vertex)
        hc = hash(str(self.crs))
        return ht + (hd << 1) + hd + hc

    def __add__(self, other):
        return Multipoint([self.vertex, other.get_vertex(self.crs)], crs=self.crs)

    @property
    def __geo_interface__(self):
        p = self.properties.copy()
        p["_karta_proj4"] = self.crs.get_proj4()
        return {"type": "Feature",
                "geometry": self.geomdict,
                "properties": p}

    @property
    def geomdict(self):
        return {"type" : "Point", "coordinates" : self.vertex}

    @property
    def x(self):
        return self.vertex[0]

    @property
    def y(self):
        return self.vertex[1]

    @property
    def z(self):
        return self.vertex[2]

    def get_vertex(self, crs=None):
        """ Return the Point vertex as a tuple. """
        if crs is None or crs==self.crs:
            return self.vertex
        else:
            return _reproject((self.x, self.y), self.crs, crs)

    def coordsxy(self, convert_to=False):
        """ Returns the x,y coordinates. Convert_to may be set to 'deg'
        or 'rad' for convenience.  """
        if convert_to == 'rad':
            return (self.x*3.14159/180., self.y*3.14159/180.)
        elif convert_to == 'deg':
            return (self.x/3.14159*180., self.y/3.14159*180.)
        else:
            return (self.x, self.y)

    def azimuth(self, other, projected=True):
        """ Returns the compass azimuth from self to other in degrees, measured
        clockwise with north at 0.

        Parameters
        ----------
        other : Point
            second point defining direction
        projected : bool, optional
            If True and self.crs is a ProjectedCRS, return the azimuth in the
            projected cooridnate system. If False, return the geodetic azimuth.
            if self.crs is a GeographicalCRS, result is always geodetic and this
            option has no effect.

        Returns
        -------
        float

        Notes
        -----
        - Return value is NaN if points are coincident
        - If CRS is geographical, returns azimuth as defined by the CRS
          instance.
        """
        x0, y0 = self.x, self.y
        if self.crs != other.crs:
            x1, y1 = other.get_vertex(self.crs)[:2]
        else:
            x1, y1 = other.x, other.y

        if (x0, y0) == (x1, y1):
            az = np.nan
        elif projected and not isinstance(self.crs, GeographicalCRS):
            az = 90.0 - math.atan2(y1-y0, x1-x0)*180.0/math.pi
            az = (az+180) % 360 - 180
        else:
            lon0, lat0 = self.crs.project(x0, y0, inverse=True)
            lon1, lat1 = self.crs.project(x1, y1, inverse=True)
            az, _, _ = self.crs.inverse(lon0, lat0, lon1, lat1)
        return az

    def walk(self, distance, azimuth, projected=True):
        """ Returns the point reached when moving in a given direction for
        a given distance from a specified starting location.

        Parameters
        ----------
        distance : float
            distance to shift point
        azimuth: float
            shift azimuth (clockwise with north at 0)
        projected: bool, optional
            If True and self.crs is a ProjectedCRS, return the compute new
            position using the projected cooridnate system. If False, return
            the geodetically correct new position. if self.crs is a
            GeographicalCRS, result is always geodetic and this option has no
            effect.
        """
        if projected and not isinstance(self.crs, GeographicalCRS):
            maz = -(azimuth+90)*math.pi/180
            maz = (90-azimuth) * math.pi/180
            x = self.x + distance*math.cos(maz)
            y = self.y + distance*math.sin(maz)
        else:
            x1, y1 = self.crs.project(self.x, self.y, inverse=True)
            x2, y2, _ = self.crs.forward(x1, y1, azimuth, distance, radians=False)
            x, y = self.crs.project(x2, y2)
        return Point((x, y), properties=self.properties, crs=self.crs)

    def distance(self, other, projected=True):
        """ Returns a distance to another Point. Distances is computed in one
        of three ways:

        1. If *self* is in a geographical coordinate system, *other* is inverse
        projected to geographical coordinates if necessary, and the geodetic
        distance is computed on the ellipsoid of *self*.

        2. If *projected* is `True`, *other* is projected to the same
        coordinate system as *self* and flat-earth distance is computed.

        3. If *projected is `False`, both *self* and *other* are inverse
        projected to geographical coordinates, and the geodetic distance is
        computed on the ellipsoid of *self*.

        If the coordinate system is geographical and a third (z) coordinate
        exists, it is assumed to have the same units as the real-world
        horizontal distance (i.e. meters).

        Parameters
        ----------
        other : Point
            point to compute distance to
        projected : bool, optional
            If True and self.crs is a ProjectedCRS, return the flat distance in
            the projected coordinate system. If False, return the ellipsoidal
            distance. if self.crs is a GeographicalCRS, result is always
            ellipsoidal/geodetic and this switch is ignored.

        Returns
        -------
        float

        Notes
        -----
        - If CRS is Geographical, returns distance as computed by the CRS
          instance (e.g. ellipsoidal or spherical).
        """
        if isinstance(self.crs, GeographicalCRS):
            lon0, lat0 = self.x, self.y
            lon1, lat1 = other.crs.project(other.x, other.y, inverse=True)
            _, _, dist = self.crs.inverse(lon0, lat0, lon1, lat1)
        elif projected:
            x0, y0 = self.vertex[:2]
            x1, y1 = other.get_vertex(self.crs)[:2]
            dist = math.sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))
        else:
            lon0, lat0 = self.crs.project(self.x, self.y, inverse=True)
            lon1, lat1 = other.crs.project(other.x, other.y, inverse=True)
            _, _, dist = self.crs.inverse(lon0, lat0, lon1, lat1)

        if 2 == len(self.vertex) == len(other.vertex):
            return dist
        else:
            return math.sqrt(dist**2. + (self.z-other.z)**2.)

    def shift(self, shift_vector, inplace=False):
        """ Shift point in space.

        Parameters
        ----------
        shift_vector : iterable
            vector with length equal to Geometry rank defining the shift
        inplace : bool
            whether shift should be in place (default False)
        """
        if len(shift_vector) != len(self.vertex):
            raise GGeoError('Shift vector length must equal geometry rank.')

        if inplace:
            self.vertex = tuple([a+b for a,b in zip(self.vertex, shift_vector)])
            return self
        else:
            vertex = tuple([a+b for a,b in zip(self.vertex, shift_vector)])
            return Point(vertex, properties=self.properties, crs=self.crs)

class MultiVertexBase(Geometry):

    #__slots__ = ["vertices", "data"]

    def __init__(self, vertices, ring=False, **kwargs):
        super(MultiVertexBase, self).__init__(**kwargs)

        if hasattr(vertices, "__next__"):
            vertices = list(vertices)

        if len(vertices) == 0:
            self.vertices = CoordString([], ring=ring)
        elif not isinstance(vertices[0], Point):
            self.vertices = CoordString(_flatten(vertices), ring=ring)
        else:
            self.vertices = CoordString([point.vertex for point in vertices], ring=ring)
            if "crs" not in kwargs:
                self.crs = vertices[0].crs

        if self.vertices.rank not in (2, 3):
            raise TypeError("Coordinates must be a sequence of sequences")
        return

    def __eq__(self, other):
        try:
            return (self._geotype == other._geotype) and \
                   (len(self.vertices) == len(other.vertices)) and \
                   np.all(np.equal(self.vertices, other.vertices)) and \
                   (self.properties == other.properties) and \
                   (self.crs == other.crs)
        except (AttributeError, TypeError, ValueError):
            return False

    def __neq__(self, other):
        return ~(self == other)

    def __hash__(self):
        ht = hash(self._geotype)
        hd = hash(self.vertices.asarray().tostring())
        hc = hash(str(self.crs))
        return ht + (hd << 1) + hd + hc

    def __getitem__(self, key):
        if isinstance(key, (int, np.int64)):
            return Point(self.vertices[key], properties=self.properties, crs=self.crs)
        elif isinstance(key, slice):
            start, stop, stride = key.indices(len(self.vertices))
            verts = self.vertices.slice(start, stop, stride)
            return type(self)(verts, properties=self.properties, crs=self.crs)
        else:
            raise KeyError('index must be integer or slice object')

    def __setitem__(self, key, value):
        if not isinstance(key, int):
            raise GGeoError('Indices must be integers')

        if getattr(value, "_geotype", None) == "Point":
            if self.crs == value.crs:
                self.vertices[key] = np.asarray(value.vertex)
            else:
                raise CRSError("Point and Geometry have different coordinate systems")
        elif len(value) == self.vertices.rank:
            self.vertices[key] = np.asarray(value)
        else:
            raise ValueError("cannot insert non-Point-like value: {0}".format(repr(value)))
        return

    def __delitem__(self, key):
        if len(self) > key:
            self.vertices = CoordString([v for i,v in enumerate(self.vertices)
                                         if i != key])
        else:
            raise GGeoError('Index ({0}) exceeds length'
                            '({1})'.format(key, len(self)))
        return

    def __iter__(self):
        return (self[i] for i in range(len(self)))

    def __len__(self):
        return len(self.vertices)

class MultiVertexMixin(object):

    def get_coordinate_lists(self, crs=None):
        """ Return horizontal coordinate lists.

        Parameters
        ----------
        crs : karta.CRS, optional
            coordinate system of output vertices
        """
        x, y = self.vertices.vectors()[:2]
        if crs is not None and (crs != self.crs):
            x, y = _reproject((x,y), self.crs, crs)
        return x, y

    @property
    def coordinates(self):
        return self.get_coordinate_lists()

    def get_vertices(self, crs=None):
        """ Return vertices as an array.

        Parameters
        ----------
        crs : karta.CRS, optional
            coordinate system of output vertices
        """
        if (crs is None) or (crs is self.crs):
            return np.array(self.vertices)
        else:
            vertices = [_reproject(v[:2], self.crs, crs)
                        for v in self.vertices]
            return np.array(vertices)
    @property
    def bbox(self):
        return self.get_bbox()

    @cache_decorator("bbox")
    def get_bbox(self, crs=None):
        """ Return the bounding box.

        Parameters
        ----------
        crs : karta.CRS
            coordinate system of output bounding box

        Returns
        -------
        tuple
            (xmin, ymin, xmax, ymax)

        Notes
        -----
        - If CRS is Geographical, returns bbox computed using a spherical
          approximation.
        """
        if crs is not None and (crs != self.crs):
            cs = CoordString(list(zip(
                *self.crs.transform(crs, *self.vertices.vectors(drop_z=True)))))
        else:
            cs = self.vertices
            crs = self.crs

        if isinstance(crs, GeographicalCRS):
            return _cdateline.bbox(cs)
        else:
            return _cvectorgeo.bbox(cs)

    @property
    def extent(self):
        return self.get_extent()

    @cache_decorator("extent")
    def get_extent(self, crs=None):
        """ Calculate geometry extent.

        Parameters
        ----------
        crs : karta.CRS
            coordinate system of output

        Returns
        -------
        tuple
            (xmin, xmax, ymin, ymax)
        """
        bb = self.get_bbox(crs=crs)
        return bb[0], bb[2], bb[1], bb[3]

    def _bbox_overlap(self, other):
        """ Whether bounding box overlaps with that of another Geometry. """
        reg0 = self.bbox
        reg1 = other.bbox
        return (reg0[0] <= reg1[2] and reg1[0] <= reg0[2] and
                reg0[1] <= reg1[3] and reg1[1] <= reg0[3])

    def shift(self, shift_vector, inplace=False):
        """ Shift geometry in space.

        Parameters
        ----------
        shift_vector : iterable
            vector with length equal to Geometry rank defining the shift
        inplace : bool
            whether shift should be in place (default False)
        """
        if len(self.vertices) == 0:
            raise GGeoError('cannot shift zero length geometry')
        if len(shift_vector) != len(self.vertices[0]):
            raise GGeoError('shift vector length must equal geometry rank')

        if inplace:
            self.vertices = CoordString(self.vertices.asarray() + shift_vector)
            self._cache = {}
            return self
        else:
            vertices = self.vertices.asarray() + shift_vector
            kw = {"properties": self.properties, "crs": self.crs}
            if hasattr(self, "data"):
                kw["data"] = self.data
            return type(self)(vertices, **kw)

    def rotate2d(self, thetad, origin=(0, 0)):
        """ Rotate rank 2 geometry.

        Parameters
        ----------
        thetad : float
            degreesof rotation
        origin : tuple of two floats, optional
            pivot for rotation (default (0, 0))
        """
        # First, shift by the origin
        self.shift([-a for a in origin])

        # Multiply by a rotation matrix
        theta = thetad / 180.0 * math.pi
        R = ((math.cos(theta), -math.sin(theta)),
             (math.sin(theta), math.cos(theta)))
        rvertices = [_matmult(R, x) for x in self.vertices]
        self.vertices = rvertices

        # Shift back
        self.shift(origin)
        self._cache = {}
        return self

    def apply_affine_transform(self, M):
        """ Apply an affine transform given by matrix *M* to data and return a
        new geometry. """
        vertices = []
        for x,y in self.get_vertices():
            vertices.append(tuple(np.dot(M, [x, y, 1])[:2]))
        return type(self)(vertices, properties=self.properties, crs=self.crs)

    def _subset(self, idxs):
        """ Return a subset defined by indices. """
        vertices = [self.vertices[i] for i in idxs]
        if hasattr(self, "data"):
            data = Table(data=[self.data._data[i] for i in idxs], fields=self.data.fields)
            return type(self)(vertices, properties=self.properties, data=data, crs=self.crs)
        else:
            return type(self)(vertices, properties=self.properties, crs=self.crs)

    def flat_distances_to(self, pt):
        """ Return the "flat Earth" distance from each vertex to a point. """
        A = np.array(self.vertices)
        P = np.tile(np.array(pt.vertex), (A.shape[0], 1))
        d = np.sqrt(np.sum((A-P)**2, 1))
        return d

    def distances_to(self, pt):
        """ Return the distance from each vertex to a point. """
        d = [pt.distance(a) for a in self]
        return np.array(d)

    def nearest_vertex_to(self, point):
        """ Returns the index of the vertex that is nearest to a point. If two
        points are equidistant, only one will be returned.

        Parameters
        ----------
        point : Point
            target point

        Returns
        -------
        int
        """
        distances = self.distances_to(point)
        idx = np.argmin(distances)
        return idx

    def any_within_poly(self, poly):
        """ Return whether any vertices are inside *poly* """
        for pt in self:
            if poly.contains(pt):
                return True
        return False

    def convex_hull(self):
        """ Return a Polygon representing the convex hull.

        Notes
        -----
        - If CRS is Geographical, returns hull computed using a spherical
          approximation. Failure may occur if the vertices are not in euclidian
          position.
        """
        if isinstance(self.crs, GeographicalCRS):
            indices = _cconvexhull.convexhull_sph(self.vertices)
        else:
            indices = _cconvexhull.convexhull(self.vertices)
        return Polygon([self.vertices[i] for i in indices], crs=self.crs)

    def to_xyfile(self, fnm, fields=None, delimiter=' ', header=None):
        """ Write data to a delimited ASCII table.

        Parameters
        ----------
        fnm : string
            filename to write to
        fields : list
            specify the fields to be written (optional, default all)

        Additional kwargs are passed to `xyfile.write_xy`.
        """
        xyfile.write_xy(self.get_vertices(), fnm, delimiter=delimiter, header=header)
        return


class ConnectedMultiVertexMixin(MultiVertexMixin):

    @cache_decorator("bbox")
    def get_bbox(self, crs=None):
        """ Dateline-aware get_bbox for geometries consisting of connected
        vertices.

        Parameters
        ----------
        crs : karta.CRS
            coordinate system of output bounding box

        Returns
        -------
        tuple
             (xmin, ymin, xmax, ymax)
        """
        if crs is not None and (crs != self.crs):
            cs = CoordString(list(zip(
                *self.crs.transform(crs, *self.vertices.vectors(drop_z=True)))))
        else:
            cs = self.vertices
            crs = self.crs

        if isinstance(crs, GeographicalCRS):
            bbox = _cdateline.bbox(cs)
        else:
            bbox = super(ConnectedMultiVertexMixin, self).get_bbox(crs=crs)
        return bbox

    @property
    @cache_decorator("length")
    def length(self):
        """ Returns the length of the line/boundary.

        Notes
        -----
        - If CRS is Geographical, returns length computed using distance
          provided by the CRS instance.
        """
        if isinstance(self.crs, GeographicalCRS):
            lon, lat = self.get_vertices()[0][:2]
            d = 0.0
            for xy in self.get_vertices()[1:]:
                d += self.crs.inverse(lon, lat, xy[0], xy[1])[2]
                lon = xy[0]
                lat = xy[1]
            return d
        else:
            return _cvectorgeo.length(self.vertices)

    @property
    def segments(self):
        """ Returns an generator of adjacent line segments. """
        return (self._subset((i,i+1)) for i in range(len(self)-1))

    @property
    def segment_tuples(self):
        """ Returns an generator of adjacent line segments as coordinate tuples. """
        return ((self.vertices[i], self.vertices[i+1])
                for i in range(len(self.vertices)-1))

    def intersects(self, other):
        """ Return whether an intersection exists with another geometry.

        Parameters
        ----------
        other : Geometry
            another geometry with multipl connected vertices

        Notes
        -----
        - If CRS is Geographical, uses a spherical approximation.
        """
        if _cintersection.bboxes_overlap(self.bbox, other.bbox):
            if isinstance(self.crs, GeographicalCRS):
                return _cintersection.intersects_sph(self.vertices, other.vertices)
            else:
                return _cintersection.intersects(self.vertices, other.vertices)

    def intersections(self, other, keep_duplicates=False):
        """ Return the intersections with another geometry as a Multipoint.

        Parameters
        ----------
        other : Geometry
            another geometry with multipl connected vertices
        keep_duplicates : bool, optional
            whether to retain duplicate intersections [default False]

        Notes
        -----
        - If CRS is Geographical, uses a spherical approximation.
        """
        if isinstance(self.crs, CartesianCRS):
            interx = _cintersection.all_intersections(self.vertices, other.vertices)
            if not keep_duplicates:
                interx = list(set(interx))
            return Multipoint(interx, crs=self.crs)
        else:
            interx = (geodesy.intersection_spherical(a, b)
                            for a in self.segment_tuples
                            for b in other.segment_tuples)
            if not keep_duplicates:
                interx = list(set(interx))
            return Multipoint(interx, crs=self.crs)

    def _nearest_to_point(self, point):
        """ Return a tuple of the shortest distance on the geometry boundary to
        a point, and the vertex at that location.

        If necessary, project coordinates to the local coordinate system.

        Parameters
        ----------
        point : Point

        Notes
        -----
        - If CRS is Geographical, uses distance defined by the CRS instance.
        """
        ptvertex = point.get_vertex(crs=self.crs)
        segments = zip(self.vertices.slice(0, -1), self.vertices.slice(1, 0))

        if isinstance(self.crs, CartesianCRS):
            func = _cvectorgeo.pt_nearest_planar
            def func(seg):
                return _cvectorgeo.pt_nearest_planar(ptvertex[0], ptvertex[1],
                                    seg[0][0], seg[0][1], seg[1][0], seg[1][1])
        else:
            fwd = self.crs.forward
            inv = self.crs.inverse
            def func(seg):
                return _cvectorgeo.pt_nearest_proj(fwd, inv, ptvertex,
                                                   seg[0], seg[1], tol=0.01)

        point_dist = map(func, segments)
        min_point = None
        min_dist = -1.0
        for i, (point, dist) in enumerate(point_dist):
            if dist < min_dist or (i == 0):
                min_point = point
                min_dist = dist

        return min_dist, min_point

    def shortest_distance_to(self, pt):
        """ Return the shortest distance from any position on the geometry
        boundary to a point.

        Parameters
        ----------
        point : Point

        Notes
        -----
        - If CRS is Geographical, uses distance defined by the CRS instance.
        """
        return self._nearest_to_point(pt)[0]

    def nearest_on_boundary(self, point):
        """ Returns the position on the geometry boundary that is nearest to
        a point. If two points are equidistant, only one will be returned.

        Parameters
        ----------
        point : Point

        Notes
        -----
        - If CRS is Geographical, uses distance defined by the CRS instance.
        """
        _, minpt = self._nearest_to_point(point)
        return Point(minpt, crs=self.crs)

    def within_distance(self, point, distance):
        """ Test whether a point is within *distance* geometry.

        Parameters
        ----------
        point : Point
        distance : float

        Notes
        -----
        - If CRS is Geographical, uses distance defined by the CRS instance.
        """
        return all(distance >= seg.shortest_distance_to(point)
                    for seg in self.segments)

    def crosses_dateline(self):
        """ Return a boolean that indicates whether any segment crosses the
        dateline.
        """
        if not isinstance(self.crs, GeographicalCRS):
            raise CRSError("Dateline detection only defined for geographical "
                           "coordinates")

        def _seg_crosses_dateline(seg):
            a, b = seg[0], seg[1]
            return (_sign(a.x) != _sign(b.x)) and (abs(a.x-b.x) > 180.0)

        return any(_seg_crosses_dateline(seg) for seg in self.segments)

class Line(MultiVertexBase, ConnectedMultiVertexMixin, GeoJSONOutMixin, ShapefileOutMixin):
    """ Line composed of connected vertices.

    Parameters
    ----------
    coords : list of 2-tuples or 3-tuples
        vertex coordinates
    properties : dict, optional
        geometry specific metadata
    crs : karta.CRS, optional
        (default Cartesian)
    """
    #__slots__ = []

    def __init__(self, vertices, **kwargs):
        """ Partial init function that creates a metadata attribute.
        """
        super(Line, self).__init__(vertices, **kwargs)
        self._geotype = "Line"
        return

    def __add__(self, other):
        return Multiline([self.vertices, other.get_vertices(self.crs)],
                         crs=self.crs)

    @property
    def __geo_interface__(self):
        p = self.properties.copy()
        p["_karta_proj4"] = self.crs.get_proj4()
        return {"type": "Feature",
                "geometry": self.geomdict,
                "properties": p}

    @property
    def geomdict(self):
        return {"type" : "LineString",
                "bbox" : self.bbox,
                "coordinates" : _as_nested_lists(self.vertices)}

    def extend(self, other):
        """ Combine two lines, provided that that the data formats are similar.
        """
        if len(self.vertices[0]) != len(other.vertices[0]):
            raise ValueError("Rank mismatch ({0} != "
                    "{1})".format(self.vertices.shape[1],
                                  other.vertices.shape[1]))
        if self._geotype != other._geotype:
            raise TypeError("Geometry mismatch ({0} != "
                    "{1})".format(self._geotype, other._geotype))

        self.vertices = np.vstack([self.vertices, other.vertices])
        self._cache = {}
        return self

    def cumulength(self):
        """ Returns the cumulative length by vertex.

        Notes
        -----
        - If CRS is Geographical, uses distance defined by the CRS instance.
        """
        d = [0.0]
        pta = self[0]
        for ptb in self[1:]:
            d_ = pta.distance(ptb)
            d.append(d_ + d[-1])
            pta = ptb
        return d

    def to_points(self, dx):
        """ Return equally spaced Point instances along line.

        Parameters
        ----------
        dx : float
            spacing of points
        """
        remainder = 0
        pt0 = self[0]
        vertices = [pt0.get_vertex()]

        for seg in self.segments:
            pos = 0
            az = seg[0].azimuth(seg[1])

            while pos < seg.length:
                distance_to_endpt = pt0.distance(seg[1])
                if distance_to_endpt >= dx:
                    pt1 = pt0.walk(dx - remainder, az)
                    pos += dx - remainder
                    vertices.append(pt1.get_vertex())
                    remainder = 0
                    pt0 = pt1
                else:
                    remainder = distance_to_endpt
                    pos = seg.length
                    pt0 = seg[1]
        return Multipoint(vertices, crs=self.crs)

    def to_npoints(self, n):
        """ Return *n* equally spaced Point instances along line.

        Parameters
        ----------
        n : int
            number of points to return
        """
        segments = self.segments
        Ltotal = self.cumulength()[-1]
        step = Ltotal / float(n-1)
        step_remaining = step

        vertices = [self[0].get_vertex()]
        x = 0.0
        pos = self[0]
        seg = next(segments)
        seg_remaining = seg.displacement()

        while x < Ltotal-1e-8:
            direction = seg[0].azimuth(seg[1])

            if step_remaining <= seg_remaining:
                pos = pos.walk(step_remaining, direction)
                x += step_remaining
                seg_remaining -= step_remaining
                step_remaining = step
                vertices.append(pos.get_vertex())
                seg.vertices[0] = np.array(pos.vertex, dtype=np.float64)

            else:
                pos = seg[1]
                x += seg_remaining
                step_remaining -= seg_remaining

                seg = next(segments, seg)
                seg_remaining = seg.displacement()
                # except StopIteration as e:
                #     if abs(Ltotal-x) > 1e-8:       # tolerance for endpoint
                #         raise e

        if len(vertices) == n-1:
            vertices.append(seg[-1].get_vertex())
        return Multipoint(vertices, crs=self.crs)

    def displacement(self):
        """ Returns the distance between the first and last vertex.

        Notes
        -----
        - If CRS is Geographical, uses distance defined by the CRS instance.
        """
        return self[0].distance(self[-1])

    def to_polygon(self):
        """ Returns a polygon. """
        return Polygon(self.vertices, properties=self.properties, crs=self.crs)

class Polygon(MultiVertexBase, ConnectedMultiVertexMixin, GeoJSONOutMixin, ShapefileOutMixin):
    """ Polygon, composed of a closed sequence of vertices.

    Parameters
    ----------
    coords : list of 2-tuples or 3-tuples
        vertex coordinates
    subs : list of Polygon instances, optional
        sub-polygons [default None]
    properties : dict, optional
        geometry specific metadata
    crs : karta.CRS, optional
        (default Cartesian)
    """
    #__slots__ = ["subs"]
    def __init__(self, vertices, subs=None, **kwargs):
        """ Partial init function that creates a metadata attribute.
        """
        super(Polygon, self).__init__(vertices, ring=True, **kwargs)
        self._geotype = "Polygon"
        if subs is not None:
            self.subs = list(subs)
        else:
            self.subs = []
        return

    def __getitem__(self, key):
        if isinstance(key, slice):
            start, stop, stride = key.indices(len(self.vertices))
            if len(self) != ((stop - start) // stride):
                return Line(self.vertices.slice(start, stop, stride),
                            properties=self.properties,
                            crs=self.crs)
        return super(Polygon, self).__getitem__(key)

    def __add__(self, other):
        return Multipolygon([[self.vertices], [other.get_vertices(self.crs)]],
                            crs=self.crs)

    @property
    def __geo_interface__(self):
        p = self.properties.copy()
        p["_karta_proj4"] = self.crs.get_proj4()
        return {"type": "Feature",
                "geometry": self.geomdict,
                "properties": p}

    @property
    def vertices_ring(self):
        # inefficient implementation
        vertices = [xy for xy in self.vertices]
        vertices.append(vertices[0])
        return CoordString(vertices)

    @property
    def geomdict(self):
        coords = [_as_nested_lists(self.vertices_ring)]
        for geom in self.subs:
            coords.append(_as_nested_lists(geom.vertices_ring))
        return {"type" : "Polygon",
                "bbox" : self.bbox,
                "coordinates" : coords}

    def _subset(self, idxs):
        """ Return a subset defined by index in *idxs*. """
        vertices = [self.vertices[i] for i in idxs]
        subset = Line(vertices, properties=self.properties, crs=self.crs)
        return subset

    def isclockwise(self):
        """ Return whether polygon winds clockwise around its interior. """
        s = sum((seg[1][0] - seg[0][0]) * (seg[1][1] + seg[0][1])
                        for seg in self.segment_tuples)
        return s > 0

    def ispolar(self, pole=None):
        """ Return True if polygon contains one pole. If the polygon contains
        neither or both poles, returns False.

        Parameters
        ----------
        pole : Point, optional
            (default point on a sphere at 0 longitude, 90 latitude)
        """

        if not isinstance(self.crs, GeographicalCRS):
            raise CRSError("ispolar defined only for geographical CRS")

        if pole is None:
            pole = Point((0, 90), crs=SphericalEarth)

        lon0 = geodesy.reduce_deg(self[-1].vertex[0])
        sum_angle = 0.0
        for vertex in self.vertices:
            lon1 = geodesy.reduce_deg(vertex[0])
            if _cdateline.crosses_dateline(lon0, lon1):
                sum_angle += 360.0 + lon1 - lon0
            else:
                sum_angle += lon1 - lon0
            lon0 = lon1

        return True if abs(sum_angle) > 1e-4 else False

    @property
    def segments(self):
        """ Returns a generator of adjacent line segments.
        """
        L = len(self.vertices)
        return itertools.chain((self._subset((i,i+1)) for i in range(len(self)-1)),
                               (self._subset((L-1,0)),))

    @property
    def segment_tuples(self):
        """ Returns a generator of adjacent line segments as coordinate
        tuples. """
        return ((self.vertices[i-1], self.vertices[i])
                for i in range(len(self.vertices)))

    @property
    def length(self):
        raise AttributeError("%s instance has no attribute 'length'" % type(self))

    @property
    def perimeter(self):
        """ Return the perimeter of the polygon. If there are sub-polygons,
        their perimeters are added recursively.

        Notes
        -----
        - If CRS is Geographical, uses distance defined by the CRS instance.
        """
        return sum(seg.length for seg in self.segments) + \
                sum([p.perimeter for p in self.subs])

    @property
    def area(self):
        """ Return the two-dimensional area of the polygon, excluding
        sub-polygons.

        Notes
        -----
        - If CRS is Geographical, uses either a spherical or an ellipsoidal
          calculation.
        """
        if isinstance(self.crs, GeographicalCRS):
            major_axis = self.crs.ellipsoid.a
            minor_axis = self.crs.ellipsoid.b

            area = 0.0
            if major_axis == minor_axis:    # Sphere
                for seg in self.segment_tuples:
                    x1, y1 = seg[0]
                    x2, y2 = seg[1]
                    area += geodesy.spherical_area(major_axis, x1, y1, x2, y2)

            else:
                for seg in self.segment_tuples:
                    x1, y1 = seg[0]
                    x2, y2 = seg[1]
                    area += geodesy.ellipsoidal_area(major_axis, minor_axis,
                                                     x1, y1, x2, y2)

        else:
            # Cartesian coordinate systems
            x, y = self.coordinates
            x0 = np.min(x)
            area = (0.5*(x[0] + x[-1]) - x0) * (y[0] - y[-1])
            area += sum((0.5*(x[i+1]+x[i]) - x0) * (y[i+1] - y[i]) for i in range(len(x)-1))
        return abs(area) - sum(sub.area for sub in self.subs)

    @property
    def centroid(self):
        """ Return Polygon centroid as a Point, ignoring sub-polygons. """
        x, y = self.coordinates
        A = 0.5 * sum(x[i]*y[i+1] - x[i+1]*y[i] for i in range(-1, len(self)-1))
        cx = sum((x[i] + x[i+1]) * (x[i]*y[i+1] - x[i+1]*y[i])
                    for i in range(-1, len(self)-1)) / (6*A)
        cy = sum((y[i] + y[i+1]) * (x[i]*y[i+1] - x[i+1]*y[i])
                    for i in range(-1, len(self)-1)) / (6*A)
        return Point((cx, cy), properties=self.properties, crs=self.crs)

    def contains(self, point):
        """ Returns True if point is inside or on the boundary of the polygon,
        and False otherwise. Uses a crossing number scheme.

        Notes
        -----
        - When the polygon is polar in a geographical coordinate system, a less
          efficient algorithm is used. For better performance, consider
          projecting to an appropriate coordinate system such as NSIDCNorth or
          NSIDCSouth beforehand.
        - Otherwise, a planar algorithm is used.
        """
        x, y = point.get_vertex(crs=self.crs)[:2]
        if isinstance(self.crs, GeographicalCRS) and self.ispolar():
            return _ccontains.contains_proj(x, y, self.vertices, self.crs) \
                    and not any(p.contains(point) for p in self.subs)
        else:
            return _ccontains.contains(x, y, self.vertices) \
                    and not any(p.contains(point) for p in self.subs)

    def to_line(self):
        """ Returns a self-closing polyline. Discards sub-polygons. """
        v = self.vertices + self.vertices[0]
        return Line(v, properties=self.properties, crs=self.crs)

class Multipart(Geometry):
    """ Base for objects consisting of multiple singular types. """

    def __init__(self, vertices, data=None, **kwargs):
        super(Multipart, self).__init__(**kwargs)

        if (data is None):
            self.data = Table(size=len(self.vertices))
        else:
            self.data = Table(data)
        return

    @property
    def d(self):
        return Indexer(self.data)

    def __eq__(self, other):
        try:
            return (self._geotype == other._geotype) and \
                   (len(self.vertices) == len(other.vertices)) and \
                   np.all(np.equal(self.vertices, other.vertices)) and \
                   (self.data == other.data) and \
                   (self.properties == other.properties) and \
                   (self.crs == other.crs)
        except (AttributeError, TypeError):
            return False

    def __neq__(self, other):
        return ~(self == other)

    def __hash__(self):
        ht = hash(self._geotype)
        hd = hash(self.vertices.asarray().tostring())
        hc = hash(str(self.crs))
        return ht + (hd << 1) + hd + hc

    def __contains__(self, other):
        if other in (part for part in self):
            return True
        else:
            return False

    def __len__(self):
        return len(self.vertices)

class Multipoint(Multipart, MultiVertexMixin, GeoJSONOutMixin, ShapefileOutMixin):
    """ Point cloud with associated attributes.

    Parameters
    ----------
    coords : list
        list of 2-tuples or 3-tuples defining vertices
    data : list, dict, Table object, or None
        point-specific data [default None]
    properties : dict or None
        geometry specific data [default None]
    crs : karta.crs.CRS subclass
        [default Cartesian]
    """

    #__slots__ = []

    def __init__(self, vertices, build_index=True, **kwargs):
        if hasattr(vertices, "__next__"):
            vertices = list(vertices)

        if len(vertices) == 0:
            self.vertices = CoordString([])
        elif not isinstance(vertices[0], Point):
            self.vertices = CoordString(vertices)
        else:
            self.vertices = CoordString([point.vertex for point in vertices])
            kwargs.setdefault("crs", vertices[0].crs)

        super(Multipoint, self).__init__(vertices, **kwargs)
        if build_index:
            self.quadtree = QuadTree(self.vertices)
        self._geotype = "Multipoint"
        return

    def __getitem__(self, key):
        if isinstance(key, numbers.Integral):
            p = self.d[key]
            p.update(self.properties)
            return Point(self.vertices[key], properties=p, crs=self.crs)
        elif isinstance(key, slice):
            start, stop, stride = key.indices(len(self.vertices))
            return Multipoint(self.vertices.slice(start, stop, stride),
                              properties=self.properties,
                              data=self.d[key], crs=self.crs)
        else:
            raise KeyError(type(key))

    def __setitem__(self, key, value):
        if isinstance(key, numbers.Integral):
            if hasattr(value, "vertex"):
                self.vertices[key] = np.array(value.get_vertex(self.crs), dtype=np.float64)
                row = []
                for field in self.data.fields:
                    row.append(value.properties.get(field, None))
                self.data[key] = tuple(row)
            else:
                if len(value) == len(self.vertices[0]):
                    verts = self.vertices.asarray()
                    verts[key] = value
                    self.vertices = CoordString(verts)
                self.data[key] = (None for _ in self.data.fields)

    @property
    def geomdict(self):
        return {"type" : "MultiPoint",
                "bbox" : self.bbox,
                "coordinates" : _as_nested_lists(self.vertices)}

    @property
    def __geo_interface__(self):
        p = self.properties.copy()
        p["_karta_proj4"] = self.crs.get_proj4()
        return {"type": "Feature",
                "geometry": self.geomdict,
                "properties": p}

    def within_radius(self, point, radius):
        """ Return subset of Multipoint within a radius. Items on the border
        are excluded.

        Parameters
        ----------
        point : Point
            point to to center filter at
        radius : float
            maximum distance from *point*

        Returns
        -------
        Multipoint
        """
        if hasattr(self, "quadtree"):
            search_bbox = (point.x-radius, point.y-radius,
                           point.x+radius, point.y+radius)
            candidate_indices = self.quadtree.search_within(point.x-radius,
                                                            point.y-radius,
                                                            point.x+radius,
                                                            point.y+radius)
            confirmed_indices = []
            for i in candidate_indices:
                if point.distance(self[i]) < radius:
                    confirmed_indices.append(i)
            confirmed_indices.sort()
        else:
            confirmed_indices = [i for i,d in enumerate(self.distances_to(point))
                                   if d < radius]
        return self._subset(confirmed_indices)

    def within_bbox(self, bbox):
        """ Return Multipoint subset that is within a square bounding box
        given by (xmin, xymin, xmax, ymax).
        """
        if hasattr(self, "quadtree"):
            indices = self.quadtree.search_within(*bbox)
            indices.sort()
        else:
            indices = [i for (i, pt) in enumerate(self)
                         if (bbox[0] < pt.x < bbox[2]) and (bbox[1] < pt.y < bbox[3])]
        return self._subset(indices)

    def within_polygon(self, poly):
        """ Return Multipoint subset that is within a polygon.
        """
        if hasattr(self, "quadtree"):
            bbox = poly.get_bbox(crs=self.crs)
            candidate_indices = self.quadtree.search_within(*bbox)
            confirmed_indices = []
            for i in candidate_indices:
                if poly.contains(self[i]):
                    confirmed_indices.append(i)
            confirmed_indices.sort()
        else:
            confirmed_indices = [i for (i, point) in enumerate(self)
                                   if poly.contains(point)]
        return self._subset(confirmed_indices)

class MultiVertexMultipartMixin(object):
    """ Mix-in class for multipart classes for which it is reasonable to ask
    whether member geometries are within or touching a multi-vertex geometry.

    E.g.
    - "which members touch this Line/Polygon?"
    - "which members are contained by this Polygon?"
    """
    @cache_decorator("bbox")
    def get_bbox(self, crs=None):
        bbs = [part.get_bbox(crs=crs) for part in self]
        xmin = min([bb[0] for bb in bbs])
        ymin = min([bb[1] for bb in bbs])
        xmax = max([bb[2] for bb in bbs])
        ymax = max([bb[3] for bb in bbs])
        return (xmin, ymin, xmax, ymax)

    @cache_decorator("extent")
    def get_extent(self, crs=None):
        bb = self.get_bbox(crs=crs)
        return bb[0], bb[2], bb[1], bb[3]

    @property
    def bbox(self):
        return self.get_bbox()

    @property
    def extent(self):
        return self.get_extent()

    def within_bbox(self, bbox, max_results=-1):
        """ Return Multipart geometry representing member geometries that are
        contained by a bounding box.

        Parameters
        ----------
        bbox : tuple
            (xmin, ymin, xmax, ymax)
        """
        indices = self.rtree.search_within(bbox, max_results=max_results)
        return type(self)([self[i] for i in indices])

    def touching_bbox(self, bbox, max_results=-1):
        """ Return Multipart geometry representing member geometries that touch
        a bounding box.

        Parameters
        ----------
        bbox : tuple
            (xmin, ymin, xmax, ymax)
        """
        indices = self.rtree.search_overlapping(bbox, max_results=max_results)
        return type(self)([self[i] for i in indices])

    def touching(self, geom):
        """ Return a Multipart geometry representing member geometries that
        touch a Line or Polygon.

        Touching is defined as intersecting a Line or Polygon boundary, or
        being contained within a Polygon.

        Parameters
        ----------
        geom : Line or Polygon

        Returns
        -------
        """
        indices = self.rtree.search_overlapping(geom.bbox)
        results = []
        if isinstance(geom, Line):
            for i in indices:
                test_geom = self[i]
                if geom.intersects(test_geom):
                    results.append(test_geom)
        elif isinstance(geom, Polygon):
            for i in indices:
                test_geom = self[i]
                pt = test_geom[0]
                if geom.contains(pt) or geom.intersects(test_geom):
                    results.append(test_geom)
        else:
            raise TypeError("argument must be Line or Polygon")
        return type(self)(results)

    def within(self, geom):
        """ Return a Multipart geometry representing member geometries
        contained within a Polygon.

        Parameters
        ----------
        geom : Polygon
        """
        if not isinstance(geom, Polygon):
            raise TypeError("argument must be Polygon")
        indices = self.rtree.search_overlapping(geom.bbox)
        results = []
        for i in indices:
            test_geom = self[i]
            pt = test_geom[0]
            if geom.contains(pt) and not geom.intersects(test_geom):
                results.append(test_geom)
        return type(self)(results)

class Multiline(Multipart, MultiVertexMultipartMixin, GeoJSONOutMixin,
                ShapefileOutMixin):
    """ Collection of lines with associated attributes.

    Parameters
    ----------
    coords : list
        list of lists of 2-tuples or 3-tuples defining lines
    data : list, dict, Table object, or None
        point-specific data [default None]
    properties : dict or None
        geometry specific data [default None]
    crs : karta.crs.CRS subclass
        [default Cartesian]
    """

    def __init__(self, vertices, build_index=True, **kwargs):
        if len(vertices) == 0:
            self.vertices = []
        elif isinstance(vertices[0], Line):
            self.vertices = [line.vertices for line in vertices]
            kwargs.setdefault("crs", vertices[0].crs)
        else:
            self.vertices = [CoordString(part) for part in vertices]
        super(Multiline, self).__init__(vertices, **kwargs)
        if build_index:
            self.rtree = RTree(self.vertices)
        self._geotype = "Multiline"
        return

    def __getitem__(self, key):
        if isinstance(key, numbers.Integral):
            properties = self.d[key]
            properties.update(self.properties)
            return Line(self.vertices[key], properties=properties, crs=self.crs)
        elif isinstance(key, slice):
            return Multiline(self.vertices[key], properties=self.properties,
                             data=self.d[key], crs=self.crs)
        else:
            raise KeyError(type(key))

    @property
    def geomdict(self):
        return {"type" : "MultiLineString",
                "bbox" : self.bbox,
                "coordinates" : _as_nested_lists(self.vertices)}

    @property
    def __geo_interface__(self):
        p = dict(_karta_proj4=self.crs.get_proj4())
        return {"type": "Feature",
                "geometry": self.geomdict,
                "properties": p}

    def get_vertices(self, crs=None):
        """ Return vertices as a list of arrays.

        Parameters
        ----------
        crs : karta.CRS, optional
            coordinate system of output vertices
        """
        if crs is None:
            return [np.array(v) for v in self.vertices]
        else:
            vertices = []
            for line in self.vertices:
                line_vertices = [_reproject(v[:2], self.crs, crs) for v in line]
                vertices.append(np.array(line_vertices))
            return vertices

    def get_coordinate_lists(self, crs=None):
        """ Returns a list of 2xn arrays representing lists of coordinates """
        ret = []
        for line_vertices in self.get_vertices(crs=crs):
            ret.append(line_vertices.T)
        return ret

class Multipolygon(Multipart, MultiVertexMultipartMixin, GeoJSONOutMixin,
                   ShapefileOutMixin):
    """ Collection of polygons with associated attributes.

    Parameters
    ----------
    coords : list
        list of lists of polygon rings, each consisting of 2-tuples or 3-tuples
    data : list, dict, Table object, or None
        point-specific data [default None]
    properties : dict or None
        geometry specific data [default None]
    crs : karta.crs.CRS subclass
        [default Cartesian]
    """

    def __init__(self, vertices, build_index=True, **kwargs):
        if len(vertices) == 0:
            self.vertices = []
        elif isinstance(vertices[0], Polygon):
            self.vertices = []
            for polygon in vertices:
                rings = [polygon.vertices]
                for sub in polygon.subs:
                    rings.append(sub)
                self.vertices.append(rings)
            kwargs.setdefault("crs", vertices[0].crs)
        else:
            self.vertices = []
            for part in vertices:
                rings = [CoordString(ring) for ring in part]
                self.vertices.append(rings)
        super(Multipolygon, self).__init__(vertices, **kwargs)
        if build_index:
            self.rtree = RTree([v[0] for v in self.vertices])
        self._geotype = "Multipolygon"
        return

    def __getitem__(self, key):
        if isinstance(key, numbers.Integral):
            properties = self.d[key]
            properties.update(self.properties)
            subs = []
            for vertices in self.vertices[key][1:]:
                subs.append(Polygon(vertices, properties=properties, crs=self.crs))
            vertices = self.vertices[key][0]
            return Polygon(vertices, subs=subs, properties=properties, crs=self.crs)
        elif isinstance(key, slice):
            return Multipolygon(self.vertices[key], properties=self.properties,
                                data=self.d[key], crs=self.crs)
        else:
            raise KeyError(type(key))

    @property
    def geomdict(self):
        return {"type" : "MultiPolygon",
                "bbox" : self.bbox,
                "coordinates" : _as_nested_lists(self.vertices)}

    @property
    def __geo_interface__(self):
        p = dict(_karta_proj4=self.crs.get_proj4())
        return {"type": "Feature",
                "geometry": self.geomdict,
                "properties": p}

    @property
    def vertices_ring(self):
        # inefficient implementation
        vertices = []
        for poly in self.vertices:
            rings = []
            for ring in poly:
                xys = [xy for xy in ring]
                xys.append(ring[0])
                rings.append(CoordString(xys))
            vertices.append(rings)
        return vertices

    def get_vertices(self, crs=None):
        """ Return vertices as a list of arrays.

        Parameters
        ----------
        crs : karta.CRS, optional
            coordinate system of output vertices
        """
        if crs is None:
            vertices = []
            for poly_vertices in self.vertices:
                vertices.append([np.array(v) for v in poly_vertices])
            return vertices
        else:
            vertices = []
            for poly_vertices in self.vertices:
                poly = []
                for ring_vertices in poly_vertices:
                    poly.append(np.array([_reproject(v[:2], self.crs, crs)
                                          for v in ring_vertices]))
                vertices.append(poly)
            return vertices

    def get_coordinate_lists(self, crs=None):
        """ Returns a list of 2xn arrays representing lists of coordinates """
        ret = []
        for poly_vertices in self.get_vertices(crs=crs):
            poly = []
            for ring_vertices in poly_vertices:
                poly.append(ring_vertices.T)
            ret.append(poly)
        return ret

def _signcross(a, b):
    """ Return sign of 2D cross product a x b """
    c = (a[0]*b[1]) - (a[1]*b[0])
    if c != 0:
        return c/abs(c)
    else:
        return 0

def _matmult(A, x):
    """ Return product of matrix A and vector x """
    b = []
    for a in A:
        b.append(sum([ai * xi for ai, xi in zip(a, x)]))
    return b

def _sign(a):
    """ Return the sign of *a* """
    if a == 0.0:
        return 1
    else:
        return a/abs(a)

def multipart_from_singleparts(parts, crs=None):
    """ Merge singlepart geometries into a multipart geometry.
    Properties contained by all inputs are stored in Multipart data attribute.

    Parameters
    ----------
    parts : iterable of singlepart Geometry instances
        e.g. list of Points, Lines, or Polygons
    crs : karta.CRS
        coordinate system of output Geometry

    Returns
    -------
    Multipart
        e.g. Multipoint, Multiline, or Multipolygon
    """
    if len(parts) == 0:
        raise ValueError("cannot construct multipart from zero singleparts")

    if crs is None:
        crs = parts[0].crs

    keys = list(parts[0].properties.keys())
    for part in parts[1:]:
        for key in keys:
            if key not in part.properties:
                keys.pop(keys.index(key))

    if len(keys) == 0:
        data = None
    else:
        data = {}
        for key in keys:
            data[key] = [part.properties[key] for part in parts]

    gt = parts[0]._geotype
    if gt == "Point":
        cls = Multipoint
        vertices = [part.get_vertex(crs=crs) for part in parts]
    elif gt == "Line":
        cls = Multiline
        vertices = [part.get_vertices(crs=crs) for part in parts]
    elif gt == "Polygon":
        cls = Multipolygon
        vertices = []
        for part in parts:
            part_vertices = [part.get_vertices(crs=crs)]
            for sub in part.subs:
                part_vertices.append(sub.get_vertices(crs=crs))
        vertices.append(part_vertices)
    else:
        raise GeometryError("cannot convert type '{0}' to multipart".format(gt))

    return cls(vertices, data=data, crs=crs)

def affine_matrix(mpa, mpb):
    """ Compute the affine transformation matrix that best matches two
    Multipoint geometries using a least squares fit.

    Output is relative to the coordinate system of the first geometry, if they
    differ.

    Parameters
    ----------
    mpa, mpb : Multipoint
        matching length collection of control points to match
    """
    if len(mpa) != len(mpb):
        raise GeometryError("Input geometries must have identical length")
    vecp = np.asarray(mpb.get_vertices(mpa.crs)).ravel()
    A = np.empty([2*len(mpa), 6], dtype=np.float64)
    for i, (x, y) in enumerate(mpa.get_vertices()):
        A[2*i:2*i+2,:] = np.kron(np.eye(2), [x, y, 1])
    M, res, rank, singvals = np.linalg.lstsq(A, vecp)
    return np.vstack([np.reshape(M, [2, 3]), np.atleast_2d([0, 0, 1])])

def get_tile_tuple(point, zoom):
    """ Return the (z, x, y) locator for an OpenStreetMap tile containing a
    point.

    Parameters
    ----------
    point : Point
        geographic point to be contained within tile
    zoom : int
        non-negative zoom level (typically 0-18)
    """
    z = int(zoom)
    dlon = 256
    dlat = 256

    lon0, lat0 = point.crs.project(*point.vertex[:2], inverse=True)
    c = 128/math.pi * 2**z
    x0 = c * (lon0*math.pi/180+math.pi)
    y0 = c * (math.pi-math.log(math.tan(math.pi/4+lat0*math.pi/360)))

    x = int(x0 // dlon)
    y = int(y0 // dlat)
    return z, x, y

