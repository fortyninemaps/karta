"""
Geographical measurement and simple analysis module for Python. Provides
Point, Multipoint, Line, and Polygon classes, with methods for simple
measurements such as distance, area, and bearing.
"""
from __future__ import division

import math
import os
import sys
import itertools
import warnings
import numpy as np
from . import geojson
from . import xyfile
from . import shp
from .table import Table, Indexer
from . import quadtree
from . import _cvectorgeo
from .. import geodesy
from ..crs import Cartesian, CartesianCRS, ProjectedCRS, GeographicalCRS
from ..crs import SphericalEarth
from ..errors import GeometryError, GGeoError, GUnitError, GInitError, CRSError

def cache_decorator(key):
    """ Returns a method decorator that stores the result of a function
    invocation under `self._cache[key]`

    Does not update cache when *crs* is specified. """
    def wrapping_func(f):
        def replacement_func(self, *args, **kwargs):
            if key in self._cache:
                ret = self._cache[key]
            else:
                ret = f(self, *args, **kwargs)
                if "crs" not in kwargs:
                    self._cache[key] = (ret)
            return ret
        return replacement_func
    return wrapping_func

class Geometry(object):
    """ This is the abstract base class for all geometry types """

    __slots__ = ["_geotype", "properties", "crs", "_cache"]

    def __init__(self, crs=Cartesian):
        self._geotype = None
        self.properties = {}
        self.crs = crs
        self._cache = {}
        return

    @staticmethod
    def _distance(pos0, pos1):
        """ Generic method for calculating distance between positions that
        respects CRS """
        if pos0.crs == pos1.crs:
            (x0,y0) = pos0.x, pos0.y
            (x1,y1) = pos1.x, pos1.y
            if isinstance(pos0.crs, GeographicalCRS):
                ((x0,x1), (y0,y1)) = pos0.crs.project([x0, x1], [y0, y1],
                                                      inverse=True)
                _, _, dist = pos0.crs.inverse(x0, y0, x1, y1, radians=False)
            else:
                dist = math.sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0))
        else:
            raise CRSError("Positions must use the same CRS")
        return dist

    @property
    def d(self):
        return Indexer(self.data)

    def to_shapefile(self, fnm):
        """ Save line to a shapefile """
        if not fnm.endswith(".shp"):
            fnm = fnm + ".shp"
        shp.ogr_write(fnm, self.__geo_interface__)
        return

    def to_geojson(self, f, indent=None, **kwargs):
        """ Write data as a GeoJSON string to a file-like object `f`.

        Parameters
        ----------
        f : file-like object to recieve the GeoJSON string

        *kwargs* include:
        crs : coordinate reference system
        bbox : an optional bounding box tuple in the form (w,e,s,n)
        """
        try:
            if not hasattr(f, "write"):
                fobj = open(f, "w")
            else:
                fobj = f
            writer = geojson.GeoJSONWriter(self, crs=self.crs, **kwargs)
            writer.write_json(fobj, indent)
        finally:
            if not hasattr(f, "write"):
                fobj.close()
        return writer


class Point(Geometry):
    """ Point object instantiated with:

    Parameters
    ----------
    coords : 2-tuple or 3-tuple
    data : list, dict, Table object, or None
        point-specific data [default None]
    properties : dict or None
        geometry specific data [default None]
    crs : karta.crs.CRS subclass
        [default CARTESIAN]
    """
    __slots__ = ["vertex", "data", "rank"]

    def __init__(self, coords, data=None, properties=None, **kwargs):
        try:
            self.rank = 2 if len(coords) == 2 else 3
        except TypeError:
            raise TypeError("Point coordinates must be a sequence")
        super(Point, self).__init__(**kwargs)
        self._geotype = "Point"
        self.vertex = coords

        if data is None:
            self.data = data
        else:
            self.data = Table(data)

        if hasattr(properties, "keys"):
            self.properties = properties
        else:
            properties = {}
        return

    def __getitem__(self, idx):
        return self.vertex[idx]

    def __repr__(self):
        return 'Point(' + str(self.vertex) + ')'

    def __eq__(self, other):
        try:
            return (tuple(self.vertex) == tuple(other.vertex)) and \
                   (self.data == other.data) and \
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

    def azimuth(self, other, crs=SphericalEarth):
        """ Returns the compass azimuth from self to other in radians (i.e.
        clockwise, with north at 0). Returns NaN if points are coincident. """

        if (self.x, self.y) == (other.x, other.y):
            az = np.nan

        elif self.crs == other.crs:
            lon0, lat0 = self.crs.project(self.x, self.y, inverse=True)
            lon1, lat1 = self.crs.project(other.x, other.y, inverse=True)
            az, _, _ = self.crs.inverse(lon0, lat0, lon1, lat1)

        else:
            raise CRSError("Azimuth undefined for points in CRS {0} and "
                           "{1}".format(self.crs, other.crs))
        return az

    def walk(self, distance, direction, radians=False):
        """ Returns the point reached when moving in a given direction for
        a given distance from a specified starting location.

        Parameters
        ----------
        distance : float
            distance to walk
        direction : float
            walk azimuth (clockwise with "north" at 0)
        """
        xg1, yg1 = self.crs.project(self.x, self.y, inverse=True)
        xg2, yg2, _ = self.crs.forward(xg1, yg1, direction, distance,
                                        radians=radians)
        x, y = self.crs.project(xg2, yg2)
        return Point((x, y), properties=self.properties, data=self.data,
                     crs=self.crs)

    def distance(self, other):
        """ Returns a distance to another Point. If the coordinate system is
        geographical and a third (z) coordinate exists, it is assumed to have
        the same units as the real-world horizontal distance (i.e. meters). """
        if self.crs != other.crs:
            raise CRSError("Points must share the same coordinate system.")
        flat_dist = self._distance(self, other)
        if 2 == self.rank == other.rank:
            return flat_dist
        else:
            return math.sqrt(flat_dist**2. + (self.z-other.z)**2.)

    def shift(self, shift_vector, inplace=False):
        """ Shift point by the amount given by a vector. """
        if len(shift_vector) != self.rank:
            raise GGeoError('Shift vector length must equal geometry rank.')

        if inplace:
            self.vertex = [a+b for a,b in zip(self.vertex, shift_vector)]
            return self
        else:
            vertex = [a+b for a,b in zip(self.vertex, shift_vector)]
            return Point(tuple(vertex), data=self.data, properties=self.properties, crs=self.crs)

    def as_geojson(self, indent=2, **kwargs):
        """ Write data as a GeoJSON string to a file-like object `f`.

        Parameters
        ----------
        f : file-like object to recieve the GeoJSON string
        crs : karta.crs.CRS subclass
            coordinate reference system, optional
        bbox : 4-tuple
            a bounding box in the form (w,e,s,n), optional
        """
        writer = geojson.GeoJSONWriter(self, **kwargs)
        return writer.print_json(indent)


class MultipointBase(Geometry):
    """ Point cloud with associated attributes. This is a base class for the
    polyline and polygon classes. """

    __slots__ = ["vertices", "data", "rank", "quadtree"]

    def __init__(self, vertices, data=None, properties=None, **kwargs):
        """ Partial init function that establishes geometry rank and creates a
        metadata attribute.
        """
        super(MultipointBase, self).__init__(**kwargs)
        self._geotype = "MultipointBase"
        vertices = list(vertices)
        if len(vertices) > 0:

            def ispoint(a):
                return getattr(a, "_geotype", None) == "Point"

            if all(ispoint(a) for a in vertices):
                pts = vertices

                # Consolidate a list of points
                rank = pts[0].rank
                crs = pts[0].crs
                if len(pts) != 1 and any(crs != pt.crs for pt in pts[1:]):
                    raise CRSError("All points must share the same CRS")
                if len(pts) != 1 and any(rank != pt.rank for pt in pts[1:]):
                    raise GInitError("Input must have consistent rank")

                self.vertices = [pt.vertex for pt in pts]
                self.rank = rank
                self.crs = crs

                # Data
                if data is not None:
                    raise ValueError("data kweyword disallowed when "
                                     "constructing from points")
                keys = set.intersection(
                            *map(lambda pt: set(pt.properties.keys()), pts))
                if len(keys) != 0:
                    d = {}
                    for k in keys:
                        d[k] = [pt.properties[k] for pt in pts]
                    self.data = Table(d)
                else:
                    self.data = None

            else:
                # Construct from a list of positions (tuples)
                self.rank = len(vertices[0])

                if not 2 <= self.rank <= 3:
                    raise GInitError("Input must be doubles or triples")
                if any(self.rank != len(v) for v in vertices):
                    raise GInitError("Input must have consistent rank")

                self.vertices = [tuple(v) for v in vertices]

                if data is None:
                    self.data = None
                else:
                    self.data = Table(data)

        else:
            self.rank = None
            self.vertices = []

        if hasattr(properties, "keys"):
            self.properties = properties
        elif properties is None:
            self.properties = {}
        else:
            raise GInitError("value provided as 'properties' must be a hash")

        self.quadtree = None
        return

    def __repr__(self):
        if len(self) < 5:
            ppverts = str(self.vertices)
        else:
            ppverts = str(self.vertices[:2])[:-1] + "..." + str(self.vertices[-2:])[1:]
        return '{typ}({verts})>'.format(
                typ=str(type(self))[:-1], verts=ppverts)

    def __len__(self):
        return len(self.vertices)

    def __getitem__(self, key):
        d = None
        if isinstance(key, (int, np.int64)):
            if self.data is not None:
                d = Table([self.data[key]], fields=self.data._fields)
            return Point(self.vertices[key], data=d, properties=self.properties,
                         crs=self.crs)
        elif isinstance(key, slice):
            if self.data is not None:
                d = Table(self.data[key], fields=self.data._fields)
            return type(self)(self.vertices[key], data=d,
                              properties=self.properties, crs=self.crs)
        else:
            raise GGeoError('Index must be an integer or a slice object')

    def __setitem__(self, key, value):
        if not isinstance(key, int):
            raise GGeoError('Indices must be integers')
        try:
            self.vertices[key] = value.vertex
            if None not in (self.data, value.data):
                self.data[key] = tuple(value.data[value.data._fields.index(f)] for f in self._fields)
            elif self.data is not None:
                self.data[key] = (None for i in range(len(self.data._fields)))

        except AttributeError:
            self.vertices[key] = value
            if self.data is not None:
                self.data[key] = (None for i in range(len(self.data._fields)))

        if getattr(value, "_geotype", None) == "Point":
            self.vertices[key] = value.vertex
            if self.data is not None:
                self.data[key] = value.data._data[0]
        elif len(value) == self.rank:
            self.vertices[key] = value
        else:
            raise GGeoError('Cannot insert non-Pointlike value with '
                            'length != {0}'.format(self.rank))
        return

    def __delitem__(self, key):
        if len(self) > key:
            del self.vertices[key]
            del self.data[key]
        else:
            raise GGeoError('Index ({0}) exceeds length'
                            '({1})'.format(key, len(self)))
        return

    def __iter__(self):
        return (self[i] for i in range(len(self)))

    def __eq__(self, other):
        try:
            return (self._geotype == other._geotype) and \
                   (self.vertices == other.vertices) and \
                   (self.data == other.data) and \
                   (self.properties == other.properties) and \
                   (self.crs == other.crs)
        except (AttributeError, TypeError):
            return False

    def __neq__(self, other):
        return ~(self == other)

    def __hash__(self):
        ht = hash(self._geotype)
        hd = hash(tuple(self.vertices))
        hc = hash(str(self.crs))
        return ht + (hd << 1) + hd + hc

    def _bbox_overlap(self, other):
        """ Return whether bounding boxes between self and another geometry
        overlap.
        """
        reg0 = self.bbox
        reg1 = other.bbox
        return (reg0[0] <= reg1[2] and reg1[0] <= reg0[2] and
                reg0[1] <= reg1[3] and reg1[1] <= reg0[3])

    @property
    def bbox(self):
        return self.get_bbox()

    @cache_decorator("bbox")
    def get_bbox(self, crs=None):
        """ Return the extent of a bounding box as
            (xmin, ymin, xmax, ymax)
        """
        x, y = self.get_coordinate_lists(crs)
        bbox = (min(x), min(y), max(x), max(y))
        return bbox

    @property
    def coordinates(self):
        return self.get_coordinate_lists()

    def print_vertices(self):
        """ Prints an enumerated list of indices. """
        for i, vertex in enumerate(self.vertices):
            print("{0}\t{1}".format(i, vertex))

    def get_vertices(self, crs=None):
        """ Return vertices as an array. """
        if (crs is None) or (crs is self.crs):
            return np.array(self.vertices)
        else:
            vertices = [_reproject(v[:2], self.crs, crs) for v in self.vertices]
            return np.array(vertices)

    def get_coordinate_lists(self, crs=None):
        """ Return horizontal coordinate lists, optionally projected to *crs*.
        """
        if self.rank == 2:
            x, y = tuple(zip(*self.vertices))
        else:
            x, y, _ = tuple(zip(*self.vertices))
        if crs is not None and (crs != self.crs):
            x, y = _reproject((x,y), self.crs, crs)
        return x, y

    def shift(self, shift_vector, inplace=False):
        """ Shift feature by the amount given by a vector. """
        if len(shift_vector) != self.rank:
            raise GGeoError('Shift vector length must equal geometry rank.')

        if self.rank == 2:
            f = lambda pt: (pt[0] + shift_vector[0], pt[1] + shift_vector[1])
        elif self.rank == 3:
            f = lambda pt: (pt[0] + shift_vector[0], pt[1] + shift_vector[1],
                            pt[2] + shift_vector[2])
        if inplace:
            self.vertices = list(map(f, self.vertices))
            self.quadtree = None
            self._cache = {}
            return self
        else:
            vertices = list(map(f, self.vertices))
            return type(self)(vertices, data=self.data, properties=self.properties, crs=self.crs)

    def _matmult(self, A, x):
        """ Return Ax=b """
        b = []
        for a in A:
            b.append(sum([ai * xi for ai, xi in zip(a, x)]))
        return b

    def rotate2d(self, thetad, origin=(0, 0)):
        """ Rotate rank 2 Multipoint around *origin* counter-clockwise by
        *thetad* degrees. """
        # First, shift by the origin
        self.shift([-a for a in origin])

        # Multiply by a rotation matrix
        theta = thetad / 180.0 * math.pi
        R = ((math.cos(theta), -math.sin(theta)),
             (math.sin(theta), math.cos(theta)))
        rvertices = [self._matmult(R, x) for x in self.vertices]
        self.vertices = rvertices

        # Shift back
        self.shift(origin)
        self.quadtree = None
        self._cache = {}
        return self

    def apply_affine_transform(self, M):
        """ Apply an affine transform given by matrix *M* to data and return a
        new geometry. """
        vertices = []
        for x,y in self.get_vertices():
            vertices.append(tuple(np.dot(M, [x, y, 1])[:2]))
        return type(self)(vertices, data=self.data, properties=self.properties,
                          crs=self.crs)

    def _subset(self, idxs):
        """ Return a subset defined by index in *idxs*. """
        if len(idxs) == 0:
            raise ValueError("attempted to extract a zero-length subset")
        vertices = [self.vertices[i] for i in idxs]
        if self.data is not None:
            data = Table([self.data[i] for i in idxs], fields=self.data.fields)
        else:
            data = None
        subset = type(self)(vertices, data=data, properties=self.properties,
                            crs=self.crs)
        return subset

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

    def nearest_point_to(self, pt):
        """ Returns the internal Point that is nearest to *pt*. If two points
        are equidistant, only one will be returned.
        """
        distances = self.distances_to(pt)
        idx = np.argmin(distances)
        return self[idx]

    @property
    def extent(self):
        return self.get_extent()

    def get_extent(self, crs=None):
        """ Calculate a bounding box. """
        def gen_minmax(G):
            """ Get the min/max from a single pass through a generator. """
            (xmin, ymin) = (xmax, ymax) = next(G)
            for (x, y) in G:
                xmin = min(xmin, x)
                xmax = max(xmax, x)
                ymin = min(ymin, y)
                ymax = max(ymax, y)
            return xmin, xmax, ymin, ymax

        if (crs is None) or (crs == self.crs):
            xmin, xmax, ymin, ymax = gen_minmax(c[:2] for c in self.vertices)
        else:
            xmin, xmax, ymin, ymax = gen_minmax(_reproject(v[:2], self.crs, crs)
                                                for v in self.vertices)
        return xmin, xmax, ymin, ymax

    def any_within_poly(self, poly):
        """ Return whether any vertices are inside *poly* """
        for pt in self:
            if poly.contains(pt):
                return True
        return False

    def convex_hull(self):
        """ Return a Polygon representing the convex hull.

        Only implemented for Cartesian-derived coordinate systems.
        """
        if not isinstance(self.crs, CartesianCRS):
            raise errors.CRSError("convex_hull only implemented for "
                                  "cartesian and projected coordinate "
                                  "systems")

        points = [pt for pt in self]

        # Find the lowermost (left?) point
        pt0 = points[0]
        idx = 0
        for i, pt in enumerate(points[1:]):
            if (pt.y < pt0.y) or ((pt.y == pt0.y) and (pt.x < pt0.x)):
                pt0 = pt
                idx = i+1
        points.pop(idx)

        # Sort CCW relative to pt0, and drop all but farthest of any duplicates
        points.sort(key=lambda pt: pt0.distance(pt))
        points.sort(key=lambda pt: _cvectorgeo.polarangle(pt0.vertex, pt.vertex))
        alpha = -1
        drop = []
        for i,pt in enumerate(points):
            a = _cvectorgeo.polarangle(pt0.vertex, pt.vertex)
            if a == alpha:
                drop.append(i)
            else:
                alpha = a

        if len(drop) != 0:
            for i in drop[::-1]:
                points.pop(i)

        # initialize convex hull
        if len(points) == 2:
            return [pt0, points[0], points[1]]
        elif len(points) == 1:
            return [pt0, points[0]]
        else:

            S = [pt0, points[0], points[1]]
            for pt in points[2:]:
                while not _cvectorgeo.isleft(S[-2].vertex, S[-1].vertex, pt.vertex):
                    S.pop()
                S.append(pt)

        return Polygon(S, crs=self.crs)

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

    def as_geojson(self, indent=2, **kwargs):
        """ Print representation of internal data as a GeoJSON string.

        Parameters
        ----------
        indent : int
            spaces to intent JSON levels
        bbox : 4-tuple
            a bounding box in the form (w,e,s,n), optional
        crs : karta.crs.CRS subclass
        """
        writer = geojson.GeoJSONWriter(self, crs=self.crs, **kwargs)
        return writer.print_json(indent)


class Multipoint(MultipointBase):
    """ Point cloud with associated attributes. This is a base class for the
    polyline and polygon classes.

    Parameters
    ----------
    coords : list of 2-tuples or 3-tuples
    data : list, dict, Table object, or None
        point-specific data [default None]
    properties : dict or None
        geometry specific data [default None]
    crs : karta.crs.CRS subclass
        [default CARTESIAN]
    """

    __slots__ = []

    def __init__(self, vertices, data=None, properties=None, **kwargs):
        """ Partial init function that establishes geometry rank and creates a
        metadata attribute.
        """
        super(Multipoint, self).__init__(vertices, data=data,
                                                   properties=properties,
                                                   **kwargs)
        self._geotype = "Multipoint"
        return

    def __contains__(self, other):
        if other in (pt for pt in self):
            return True
        else:
            return False

    @property
    def __geo_interface__(self):
        p = self.properties.copy()
        p["_karta_proj4"] = self.crs.get_proj4()
        return {"type": "Feature",
                "geometry": self.geomdict,
                "properties": p}

    @property
    def geomdict(self):
        return {"type" : "MultiPoint",
                "bbox" : self.bbox,
                "coordinates" : self.vertices}

    def within_radius(self, pt, radius):
        """ Return Multipoint of subset that is within *radius* of *pt*.
        """
        if self.quadtree is None:
            distances = self.distances_to(pt)
            indices = [i for i,d in enumerate(distances) if d <= radius]
        else:
            search_bbox = (pt.x-radius, pt.y-radius, pt.x+radius, pt.y+radius)
            possible_pt_tuples = self.quadtree.getfrombbox(search_bbox)
            possible_pts = []
            for t in possible_pt_tuples:
                possible_pts.append(Point((t[0], t[1]), properties={"idx": t[2]},
                                                        crs=self.crs))
            indices = [p.properties["idx"] for p in possible_pts
                            if pt.distance(p) <= radius]

        return self._subset(indices)

    def within_bbox(self, bbox):
        """ Return Multipoint subset that is within a square bounding box
        given by (xmin, xymin, xmax, ymax).
        """
        filtbbox = lambda pt: (bbox[0] <= pt.vertex[0] <= bbox[2]) and \
                              (bbox[1] <= pt.vertex[1] <= bbox[3])
        indices = [i for (i, pt) in enumerate(self) if filtbbox(pt)]
        return self._subset(indices)

    def within_polygon(self, poly):
        """ Return Multipoint subset that is within a polygon.
        """
        if self.quadtree is None:
            indices = [i for (i, pt) in enumerate(self) if poly.contains(pt)]
        else:
            ext = poly.get_extent(crs=self.crs)
            search_bbox = (ext[0], ext[2], ext[1], ext[3])
            possible_pt_tuples = self.quadtree.getfrombbox(search_bbox)
            possible_pts = []
            for t in possible_pt_tuples:
                possible_pts.append(Point((t[0], t[1]), properties={"idx": t[2]},
                                                        crs=self.crs))
            indices = [p.properties["idx"] for p in possible_pts
                            if poly.contains(p)]
        return self._subset(indices)

    def build_quadtree(self, buff=1e-8):
        """ Construct an internal quadtree with the current geometry data.

        Parameters
        ----------
        buff : float or 4-tuple of floats
            specifies a spatial buffer to create around the current point
            bounding box, permitting the geometry to grow after the quadtree
            has been initialized. *buff* may be a scalar or a sequence of
            (left, right, bottom, top). (optional)
        """
        try:
            bf = (buff[0], buff[1], buff[2], buff[3])
        except TypeError:
            bf = (buff, buff, buff, buff)

        x0, y0, x1, y1 = self.bbox
        self.quadtree = quadtree.QuadTree((x0-bf[0], y0-bf[1], x1+bf[2], y1+bf[3]))
        for i, pt in enumerate(self):
            self.quadtree.addpt((pt.x, pt.y, i))
        return

class ConnectedMultipoint(MultipointBase):
    """ Class for Multipoints in which vertices are assumed to be connected. """

    __slots__ = []

    @property
    def bbox(self):
        return self.get_bbox()

    @cache_decorator("bbox")
    def get_bbox(self, crs=None):
        if isinstance(self.crs, GeographicalCRS):
            X, Y = self.get_coordinate_lists(crs)
            xmin = xmax = X[0]
            ymin = ymax = Y[0]
            rot = 0
            for i in range(len(X)-1):
                x0, x1 = X[i], X[i+1]
                ymin = min(ymin, Y[i+1])
                ymax = max(ymax, Y[i+1])
                seg = Line([(X[i], Y[i]), (X[i+1], Y[i+1])], crs=crs)
                if self._seg_crosses_dateline(seg):
                    if x0 > x1:     # east to west
                        rot += 360
                    else:           # west to east
                        rot -= 360

                    xmin = min(xmin, x1+rot)
                    xmax = max(xmax, x1+rot)
                else:
                    if x0 > x1:
                        xmin = min(xmin, x1)
                    else:
                        xmax = max(xmax, x1)

            xmin = (xmin+180) % 360 - 180
            xmax = (xmax+180) % 360 - 180
            bbox = (xmin, ymin, xmax, ymax)
        else:
            bbox = super(ConnectedMultipoint, self).get_bbox()
        return bbox

    @property
    @cache_decorator("length")
    def length(self):
        """ Returns the length of the line/boundary. """
        points = [Point(v, crs=self.crs) for v in self.vertices]
        distances = [a.distance(b) for a, b in zip(points[:-1], points[1:])]
        return sum(distances)

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
        """ Return whether an intersection exists with another geometry. """
        if isinstance(self.crs, CartesianCRS):
            interxbool = (np.nan in _cvectorgeo.intersection(a[0][0], a[1][0], b[0][0], b[1][0],
                                                             a[0][1], a[1][1], b[0][1], b[1][1])
                        for a in self.segments for b in other.segments)
            if self._bbox_overlap(other) and (False in interxbool):
                return True
            else:
                return False
        else:
            for a in self.segments:
                for b in other.segments:
                    try:
                        geodesy.intersection_spherical(a, b)
                    except geodesy.NoIntersection:
                        continue
                    return True
            return False

    def intersections(self, other, keep_duplicates=False):
        """ Return the intersections with another geometry as a Multipoint. """
        if isinstance(self.crs, CartesianCRS):
            interx = (_cvectorgeo.intersection(a[0][0], a[1][0], b[0][0], b[1][0],
                                              a[0][1], a[1][1], b[0][1], b[1][1])
                         for a in self.segments for b in other.segments)
            if not keep_duplicates:
                interx = set(interx)
            interx_points = []
            for vertex in interx:
                if np.nan not in vertex:
                    interx_points.append(Point(vertex, properties=self.properties,
                                               crs=self.crs))
            return Multipoint(interx_points)
        else:
            # FIXME: implement for ellipsoidal coordinate systems
            interx = (geodesy.intersection_spherical(a, b)
                         for a in self.segments for b in other.segments)
            if not keep_duplicates:
                interx = set(interx)
            interx_points = []
            for vertex in interx:
                if np.nan not in vertex:
                    interx_points.append(Point(vertex, properties=self.properties,
                                               crs=self.crs))
            return Multipoint(interx_points)


    def _nearest_to_point(self, pt):
        """ Return a tuple of the shortest distance on the geometry boundary to
        *pt*, and the vertex at that location.

        If necessary, project coordinates to the coordinate system of *self*.
        """
        ptvertex = pt.get_vertex(crs=self.crs)
        segments = zip(self.vertices[:-1], self.vertices[1:])

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
        minpt = None
        mindist = -1.0
        for i, (pt, d) in enumerate(point_dist):
            if d < mindist or (i == 0):
                minpt = pt
                mindist = d

        return mindist, minpt

    def shortest_distance_to(self, pt):
        """ Return the shortest distance from any position on the geometry
        boundary to *pt* (Point).
        """
        return self._nearest_to_point(pt)[0]

    def nearest_on_boundary(self, pt):
        """ Returns the position on the geometry boundary that is nearest to
        *pt* (Point). If two points are equidistant, only one will be returned.
        """
        _, minpt = self._nearest_to_point(pt)
        return Point(minpt, crs=self.crs)

    def within_distance(self, pt, distance):
        """ Test whether a point is within *distance* of a ConnectedMultipoint.
        """
        return all(distance >= seg.shortest_distance_to(pt) for seg in self.segments)

    @staticmethod
    def _seg_crosses_dateline(seg):
        a, b = seg[0], seg[1]
        return (sign(a.x) != sign(b.x)) and (abs(a.x-b.x) > 180.0)

    def crosses_dateline(self):
        """ Return a boolean that indicates whether any segment crosses the
        dateline """
        if not isinstance(self.crs, GeographicalCRS):
            raise CRSError("Dateline detection only defined for geographical "
                           "coordinates")

        return any(self._seg_crosses_dateline(seg) for seg in self.segments)


class Line(ConnectedMultipoint):
    """ Line composed of connected vertices.

    Parameters
    ----------
    coords : list of 2-tuples or 3-tuples
    data : list, dict, Table object, or None
        point-specific data [default None]
    properties : dict or None
        geometry specific data [default None]
    crs : karta.crs.CRS subclass
        [default CARTESIAN]
    """
    __slots__ = []

    def __init__(self, vertices, data=None, properties=None, **kwargs):
        """ Partial init function that establishes geometry rank and creates a
        metadata attribute.
        """
        super(ConnectedMultipoint, self).__init__(vertices, data=data,
                                                  properties=properties,
                                                  **kwargs)
        self._geotype = "Line"
        return

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
                "coordinates" : self.vertices}

    def extend(self, other):
        """ Combine two lines, provided that that the data formats are similar.
        """
        if self.rank != other.rank:
            raise GGeoError("Rank mismatch ({0} != {1})".format(self.rank, other.rank))
        if self._geotype != other._geotype:
            raise GGeoError("Geometry mismatch ({0} != {1})".format(self._geotype, other._geotype))

        if None not in (self.data, other.data):
            self.data._data.extend(other.data._data)
        elif self.data == other.data:
            self.data = None
        else:
            raise GGeoError('Cannot add geometries with mismatched metadata')
        i = len(self)
        self.vertices.extend(other.vertices)

        if self.quadtree is not None:
            for pt in other:
                self.quadtree.addpt((pt.x, pt.y, i))
                i += 1
        self._cache = {}
        return self

    def cumulength(self):
        """ Returns the cumulative length by vertex. """
        d = [0.0]
        pta = self[0]
        for ptb in self[1:]:
            d_ = pta.distance(ptb)
            d.append(d_ + d[-1])
            pta = ptb
        return d

    def subsection(self, n):
        """ Return *n* equally spaced Point instances along line. """
        segments = self.segments
        Ltotal = self.cumulength()[-1]
        step = Ltotal / float(n-1)
        step_remaining = step

        points = [self[0]]
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
                points.append(pos)
                seg.vertices[0] = pos.vertex

            else:
                pos = seg[1]
                x += seg_remaining
                step_remaining -= seg_remaining

                seg = next(segments, seg)
                seg_remaining = seg.displacement()
                # except StopIteration as e:
                #     if abs(Ltotal-x) > 1e-8:       # tolerance for endpoint
                #         raise e

        if len(points) == n-1:
            points.append(seg[-1])
        return points

    def displacement(self):
        """ Returns the distance between the first and last vertex. """
        return self[0].distance(self[-1])

    def to_polygon(self):
        """ Returns a polygon. """
        return Polygon(self.vertices, data=self.data, properties=self.properties,
                       crs=self.crs)


class Polygon(ConnectedMultipoint):
    """ Polygon, composed of a closed sequence of vertices.

    Parameters
    ----------
    coords : list of 2-tuples or 3-tuples
    data : list, dict, Table object, or None
        point-specific data [default None]
    properties : dict or None
        geometry specific data [default None]
    subs : list of Polygon instances or None
        sub-polygons [default None]
    crs : karta.crs.CRS subclass
        [default CARTESIAN]
    """
    __slots__ = ["subs"]

    def __init__(self, vertices, data=None, properties=None, subs=None, **kwargs):
        vertices = list(vertices)
        ConnectedMultipoint.__init__(self, vertices, data=data,
                                     properties=properties, **kwargs)
        self._geotype = "Polygon"
        self.subs = subs if subs is not None else []
        return

    def __getitem__(self, key):
        if isinstance(key, slice):
            ind = key.indices(len(self))
            if len(self) != ((ind[1] - ind[0]) // ind[2]):
                if self.data is None:
                    d = None
                else:
                    d = self.data[key]
                return Line(self.vertices[key], data=d,
                            properties=self.properties, crs=self.crs)
        return super(Polygon, self).__getitem__(key)

    @property
    def __geo_interface__(self):
        p = self.properties.copy()
        p["_karta_proj4"] = self.crs.get_proj4()
        return {"type": "Feature",
                "geometry": self.geomdict,
                "properties": p}

    @property
    def geomdict(self):
        coords = [self.vertices]
        for geom in self.subs:
            coords.append(geom.vertices)
        return {"type" : "Polygon",
                "bbox" : self.bbox,
                "coordinates" : coords}

    def _subset(self, idxs):
        """ Return a subset defined by index in *idxs*. """
        vertices = [self.vertices[i] for i in idxs]
        if self.data is None:
            data = None
        else:
            data = Table([self.data[i] for i in idxs], fields=self.data.fields)
        subset = Line(vertices, data=data, properties=self.properties,
                      crs=self.crs)
        return subset

    def isclockwise(self):
        """ Return whether polygon winds clockwise around its interior. """
        s = sum((seg[1][0] - seg[0][0]) * (seg[1][1] + seg[0][1])
                for seg in self.segments)
        return s > 0

    def ispolar(self, pole=None):
        """ Return True if polygon contains one pole. If the polygon contains
        neither or both poles, returns False.
        """

        if not isinstance(self.crs, GeographicalCRS):
            raise CRSError("ispolar defined only for geographical CRS")

        if pole is None:
            pole = Point((0, 90), crs=SphericalEarth)

        lon0 = geodesy.reduce_deg(self[-1].vertex[0])
        sum_angle = 0.0

        for vertex in self.vertices:

            lon1 = geodesy.reduce_deg(vertex[0])

            if (sign(lon0) == -sign(lon1)) and \
                ((min(abs(lon0-180), abs(lon0+180)) + min(abs(lon1-180), abs(lon1+180))) < (abs(lon0) + abs(lon1))):
                # Longitudes span the dateline
                sum_angle += 360.0 + lon1 - lon0
            else:
                sum_angle += lon1 - lon0

            lon0 = lon1

        return True if abs(sum_angle) > 1e-4 else False

    @property
    def segments(self):
        """ Returns an generator of adjacent line segments.
        Unique to Polygon: appends a final segment to close the Polygon.
        """
        L = len(self.vertices)
        return itertools.chain((self._subset((i,i+1)) for i in range(len(self)-1)),
                               (self._subset((L-1,0)),))

    @property
    def segment_tuples(self):
        """ Returns an generator of adjacent line segments as coordinate tuples. """
        return ((self.vertices[i-1], self.vertices[i]) for i in range(len(self.vertices)))

    @property
    def length(self):
        raise AttributeError("%s instance has no attribute 'length'" % type(self))

    @property
    def perimeter(self):
        """ Return the perimeter of the polygon. If there are sub-polygons,
        their perimeters are added recursively. """
        return sum(seg.length for seg in self.segments) + \
                sum([p.perimeter for p in self.subs])

    @property
    def area(self):
        """ Return the two-dimensional area of the polygon, excluding
        sub-polygons. """

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

    @staticmethod
    def _signcross(a, b):
        """ Return sign of 2D cross product a x b """
        c = (a[0]*b[1]) - (a[1]*b[0])
        if c != 0:
            return c/abs(c)
        else:
            return 0

    def contains(self, pt):
        """ Returns True if pt is inside or on the boundary of the polygon, and
        False otherwise. Uses a crossing number scheme.

        Behaviour may not be defined for polar geographical polygons.
        """
        if isinstance(self.crs, GeographicalCRS):
            if self.ispolar():
                raise NotImplementedError(
                    "Membership tests on polar geographical polygons not "
                    "implemented. As a workaround, transform to an appropriate "
                    "ProjectedCRS first.")

        x, y = pt.get_vertex(crs=self.crs)[:2]
        cnt = 0
        for seg in self.segment_tuples:
            (a, b) = seg
            if _cvectorgeo.intersects_cn(x, y, a[0], b[0], a[1], b[1]):
                cnt += 1
        return cnt % 2 == 1 and not any(p.contains(pt) for p in self.subs)

    def to_line(self):
        """ Returns a self-closing polyline. Discards sub-polygons. """
        v = self.vertices + self.vertices[0]
        return Line(v, properties=self.properties, data=self.data, crs=self.crs)

def _reproject(xy, crs1, crs2):
    """ Reproject a coordinate (or 2-tuple of x and y vectors) from *crs1* to
    *crs2*. """
    return crs1.transform(crs2, *xy)

def points_to_multipoint(points):
    """ Merge *points* into a Multipoint instance. Point properties are stored
    as Multipoint data. All points must use the same CRS.
    """
    crs = points[0].crs
    if not all(pt.crs == crs for pt in points):
        raise CRSError("All points must share the same CRS")

    keys = list(points[0].properties.keys())
    for pt in points[1:]:
        for key in keys:
            if key not in pt.properties:
                keys.pop(keys.index(key))

    ptdata = {}
    for key in keys:
        ptdata[key] = [pt.properties[key] for pt in points]

    vertices = [pt.vertex for pt in points]

    return Multipoint(vertices, data=ptdata, crs=crs)

def affine_matrix(mpa, mpb):
    """ Compute the affine transformation matrix that projects Multipoint mpa
    to Multipoint mpb using a least squares fit. """
    if len(mpa) != len(mpb):
        raise GeometryError("Input geometries must have identical length")
    vecp = np.asarray(mpb.get_vertices()).ravel()
    A = np.empty([2*len(mpa), 6], dtype=np.float64)
    for i, (x, y) in enumerate(mpa.get_vertices()):
        A[2*i:2*i+2,:] = np.kron(np.eye(2), [x, y, 1])
    M, res, rank, singvals = np.linalg.lstsq(A, vecp)
    return np.vstack([np.reshape(M, [2, 3]), np.atleast_2d([0, 0, 1])])

def sign(a):
    if a == 0.0:
        return 1
    else:
        return a/abs(a)

def get_tile_tuple(pt, zoom):
    """ Return the (z, x, y) locator for an OpenStreetMap tile containing *pt*.

    Parameters
    ----------
    pt : geographic point to be contained within tile (Point)
    zoom : non-negative zoom level (typically 0-18) (int)
    """
    z = int(zoom)
    ntiles = 2**z
    dlon = 256
    dlat = 256

    lon0, lat0 = pt.crs.project(*pt.vertex, inverse=True)
    c = 128/math.pi * 2**z
    x0 = c * (lon0*math.pi/180+math.pi)
    y0 = c * (pi-math.log(math.tan(math.pi/4+lat0*math.pi/360)))

    x = int(x0 // dlon)
    y = int(y0 // dlat)
    return z, x, y

