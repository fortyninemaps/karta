"""
Geographical measurement and simple analysis module for Python. Provides
Point, Multipoint, Line, and Polygon classes, with methods for simple
measurements such as distance, area, and bearing.
"""
from __future__ import division

import math
import sys
import itertools
import traceback
import numpy as np
from . import vtk
from . import geojson
from . import xyfile
from . import shp
from ..crs import Cartesian, SphericalEarth, CRSError
from .metadata import Metadata
from . import _vectorgeo

try:
    from . import _cvectorgeo as _vecgeo
except ImportError:
    sys.stderr.write("falling back on slow _vectorgeo")
    _vecgeo = _vectorgeo

class Geometry(object):
    """ This is the abstract base class for all geometry types """
    _geotype = None

    def __init__(self, crs=Cartesian):
        self.properties = {}
        self._crs = crs
        return

    @staticmethod
    def _distance(pos0, pos1):
        """ Generic method for calculating distance between positions that
        respects CRS """
        if pos0._crs == pos1._crs:
            (x0,y0) = pos0.x, pos0.y
            (x1,y1) = pos1.x, pos1.y
            ((x0,x1), (y0,y1)) = pos0._crs.project([x0, x1], [y0, y1],
                                                   inverse=True)
            _, _, dist = pos0._crs.inverse(x0, y0, x1, y1, radians=False)
        else:
            raise CRSError("Positions must use the same CRS")
        return dist

    @property
    def crs(self):
        return self._crs

    def add_property(self, name, value):
        """ Insert a property (name -> value) into the properties dict, raising
        a NameError if the property already exists. Compared to directly
        manipulating the properties dict, this avoids accidents. """
        if name not in self.properties:
            self.properties[name] = value
        else:
            raise NameError("property '{0}' already exists in {1}".format(name, type(self)))
        return

    def delete_property(self, name):
        """ Equivalent to

            del self.properties[name]
        """
        if name in self.properties:
            del self.properties[name]
        return


class Point(Geometry):
    """ Point object instantiated with:

    *coords*        2-tuple or 3-tuple
    *data*          List, dictionary, or Metadata object for point-specific
                    data [default None]
    *properties*    Dictionary of geometry specific data [default None]
    *crs*           Coordinate reference system instance [default CARTESIAN]
    """
    _geotype = "Point"

    def __init__(self, coords, data=None, properties=None, copy_metadata=True,
                 **kwargs):
        if not hasattr(coords, "__iter__"):
            raise ValueError("Point coordinates must be a sequence")
        super(Point, self).__init__(**kwargs)
        self.vertex = coords
        self.rank = 2 if len(coords) == 2 else 3

        self.data = Metadata(data, singleton=True, copydata=copy_metadata)
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
                   (self.properties == other.properties)
        except AttributeError:
            return False

    @property
    def __geo_interface__(self):
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
        if crs is None or crs==self._crs:
            return self.vertex
        else:
            vg = self.crs.project(self.x, self.y, inverse=True)
            return crs.project(*vg)

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

        if self.coordsxy() == other.coordsxy():
            az = np.nan

        elif self._crs == other._crs:
            lon0, lat0 = self._crs.project(self.x, self.y, inverse=True)
            lon1, lat1 = self._crs.project(other.x, other.y, inverse=True)
            az, _, _ = self._crs.inverse(lon0, lat0, lon1, lat1)

        else:
            raise CRSError("Azimuth undefined for points in CRS {0} and "
                           "{1}".format(self._crs, other._crs))
        return az

    def walk(self, distance, direction, radians=False):
        """ Returns the point reached when moving in a given direction for
        a given distance from a specified starting location.

            distance (float): distance to walk
            direction (float): walk azimuth (clockwise with "north" at 0)
        """
        xg1, yg1 = self._crs.project(self.x, self.y, inverse=True)
        xg2, yg2, _ = self._crs.forward(xg1, yg1, direction, distance,
                                        radians=radians)
        x, y = self.crs.project(xg2, yg2)
        return Point((x, y), properties=self.properties, data=self.data,
                     crs=self._crs, copy_metadata=False)

    def distance(self, other):
        """ Returns a distance to another Point. If the coordinate system is
        geographical and a third (z) coordinate exists, it is assumed to have
        the same units as the real-world horizontal distance (i.e. meters). """
        if self._crs != other._crs:
            raise CRSError("Points must share the same coordinate system.")
        flat_dist = self._distance(self, other)
        if 2 == self.rank == other.rank:
            return flat_dist
        else:
            return math.sqrt(flat_dist**2. + (self.z-other.z)**2.)

    def shift(self, shift_vector):
        """ Shift point by the amount given by a vector. Operation occurs
        in-place """
        if len(shift_vector) != self.rank:
            raise GGeoError('Shift vector length must equal geometry rank.')

        self.vertex = tuple([a+b for a,b in zip(self.vertex, shift_vector)])
        return self

    def as_geojson(self, **kwargs):
        """ Write data as a GeoJSON string to a file-like object `f`.

        Parameters
        ----------
        f : file-like object to recieve the GeoJSON string

        *kwargs* include:
        crs : coordinate reference system
        bbox : an optional bounding box tuple in the form (w,e,s,n)
        """
        writer = geojson.GeoJSONWriter(self, **kwargs)
        return writer.print_json()

    def to_geojson(self, f, **kwargs):
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
            writer = geojson.GeoJSONWriter(self, **kwargs)
            writer.write_json(fobj)
        finally:
            if not hasattr(f, "write"):
                fobj.close()
        return writer


class MultipointBase(Geometry):
    """ Point cloud with associated attributes. This is a base class for the
    polyline and polygon classes. """
    _geotype = "MultipointBase"

    def __init__(self, vertices, data=None, properties=None, copy_metadata=True,
                 **kwargs):
        """ Partial init function that establishes geometry rank and creates a
        metadata attribute.
        """
        super(MultipointBase, self).__init__(**kwargs)
        vertices = list(vertices)
        if len(vertices) > 0:

            def ispoint(a):
                return getattr(a, "_geotype", None) == "Point"

            if all(ispoint(a) for a in vertices):

                # Consolidate a list of points
                rank = vertices[0].rank
                crs = vertices[0].crs
                if len(vertices) != 1 and any(crs != pt._crs for pt in vertices[1:]):
                    raise CRSError("All points must share the same CRS")
                if len(vertices) != 1 and any(rank != pt.rank for pt in vertices[1:]):
                    raise GInitError("Input must have consistent rank")

                keys = list(vertices[0].data.keys())
                if len(vertices) != 1:
                    for pt in vertices[1:]:
                        for key in keys:
                            if key not in pt.data:
                                keys.pop(keys.index(key))

                ptdata = {}
                for key in keys:
                    ptdata[key] = [pt.data[key] for pt in vertices]

                if data is None:
                    data = ptdata
                else:
                    data.update(ptdata)

                self.vertices = [pt.vertex for pt in vertices]
                self.rank = rank
                self._crs = crs

            else:

                # Construct from a list of positions (tuples)

                self.rank = len(vertices[0])

                if not 2 <= self.rank <= 3:
                    raise GInitError("Input must be doubles or triples")
                if any(self.rank != len(v) for v in vertices):
                    raise GInitError("Input must have consistent rank")

                self.vertices = [tuple(v) for v in vertices]

        else:
            self.rank = None
            self.vertices = []

        if hasattr(properties, "keys"):
            self.properties = properties
        elif properties is None:
            self.properties = {}
        else:
            raise GInitError("value provided as 'properties' must be a hash")

        self.data = Metadata(data, copydata=copy_metadata)
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
        if isinstance(key, (int, np.int64)):
            return Point(self.vertices[key], data=self.data[key],
                         properties=self.properties, crs=self._crs,
                         copy_metadata=False)
        elif isinstance(key, slice):
            return type(self)(self.vertices[key], data=self.data[key],
                              properties=self.properties, crs=self._crs,
                              copy_metadata=False)
        else:
            raise GGeoError('Index must be an integer or a slice object')

    def __setitem__(self, key, value):
        if not isinstance(key, int):
            raise GGeoError('Indices must be integers')
        if getattr(value, "_geotype", None) == "Point":
            self.vertices[key] = value.vertex
            self.data[key] = list(value.data.values())[0]
        elif len(value) == self.rank:
            self.vertices[key] = value
            self.data[key] = None
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
        return hasattr(other, "_geotype") and \
               (self._geotype == other._geotype) and \
               (len(self) == len(other)) and \
               all(a==b for a,b in zip(self, other))

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
        """ Return the extents of a bounding box as
            (xmin, ymin, xmax, ymax)
        """
        x, y = self.get_coordinate_lists()
        bbox = (min(x), min(y), max(x), max(y))
        return bbox

    @property
    def coordinates(self):
        return self.get_coordinate_lists()

    def print_vertices(self):
        """ Prints an enumerated list of indices. """
        for i, vertex in enumerate(self.vertices):
            print("{0}\t{1}".format(i, vertex))

    def get_vertices(self):
        """ Return vertices as a list of tuples. """
        return np.array(self.vertices)

    def get_coordinate_lists(self, crs=None):
        """ Return horizontal coordinate lists, optionally projected to *crs*.
        """
        if self.rank == 2:
            x, y = tuple(zip(*self.vertices))
        else:
            x, y, _ = tuple(zip(*self.vertices))
        if crs is not None and (crs != self._crs):
            xg, yg = self.crs.project(x, y, inverse=True)
            x, y = crs.project(xg, yg)
        return x, y

    def shift(self, shift_vector):
        """ Shift feature by the amount given by a vector. Operation
        occurs in-place """
        if len(shift_vector) != self.rank:
            raise GGeoError('Shift vector length must equal geometry rank.')

        if self.rank == 2:
            f = lambda pt: (pt[0] + shift_vector[0], pt[1] + shift_vector[1])
        elif self.rank == 3:
            f = lambda pt: (pt[0] + shift_vector[0], pt[1] + shift_vector[1],
                            pt[2] + shift_vector[2])
        self.vertices = list(map(f, self.vertices))
        return self

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
        return self

    def apply_affine_transform(self, M):
        """ Apply an affine transform given by matrix *M* to data and return a
        new geometry. """
        vertices = []
        for x,y in self.get_vertices():
            vertices.append(tuple(np.dot(M, [x, y, 1])[:2]))
        return type(self)(vertices, data=self.data, properties=self.properties,
                          crs=self._crs)

    def _subset(self, idxs):
        """ Return a subset defined by index in *idxs*. """
        vertices = [self.vertices[i] for i in idxs]
        data = self.data.sub(idxs)
        subset = type(self)(vertices, data=data, properties=self.properties,
                            crs=self._crs, copy_metadata=False)
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

    def get_extents(self):
        """ Calculate a bounding box. """
        def gen_minmax(G):
            """ Get the min/max from a single pass through a generator. """
            mn = mx = next(G)
            for x in G:
                mn = min(mn, x)
                mx = max(mx, x)
            return mn, mx
        # Get the min/max for a generator defined for each dimension
        return list(map(gen_minmax,
                    map(lambda i: (c[i] for c in self.vertices),
                        range(self.rank))))

    def any_within_poly(self, poly):
        """ Return whether any vertices are inside *poly* """
        for pt in self:
            if poly.contains(pt):
                return True
        return False

    def convex_hull(self):
        """ Return a Polygon representing the convex hull. Assumes that the CRS
        may be treated as cartesian. """
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
        points.sort(key=lambda pt: _vecgeo.polarangle(pt0.vertex, pt.vertex))
        alpha = -1
        drop = []
        for i,pt in enumerate(points):
            a = _vecgeo.polarangle(pt0.vertex, pt.vertex)
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
                while not _vecgeo.isleft(S[-2].vertex, S[-1].vertex, pt.vertex):
                    S.pop()
                S.append(pt)

        return Polygon(S)

    def to_xyfile(self, fnm, fields=None, delimiter=' ', header=None):
        """ Write data to a delimited ASCII table.

        fnm         :   filename to write to

        kwargs:
        fields      :   specify the fields to be written (default all)

        Additional kwargs are passed to `xyfile.write_xy`.
        """
        xyfile.write_xy(self.get_vertices(), fnm, delimiter=delimiter, header=header)
        return

    def as_geojson(self, **kwargs):
        """ Print representation of internal data as a GeoJSON string.

        Parameters
        ----------
        crs : coordinate reference system
        bbox : an optional bounding box tuple in the form (w,e,s,n)
        """
        writer = geojson.GeoJSONWriter(self, crs=self._crs, **kwargs)
        return writer.print_json()

    def to_geojson(self, f, **kwargs):
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
            writer = geojson.GeoJSONWriter(self, crs=self._crs, **kwargs)
            writer.write_json(fobj)
        finally:
            if not hasattr(f, "write"):
                fobj.close()
        return writer

    def to_vtk(self, f, **kwargs):
        """ Write data to an ASCII VTK .vtp file. """
        if not hasattr(f, "write"):
            with open(f, "w") as fobj:
                vtk.mp2vtp(self, fobj, **kwargs)
        else:
            vtk.mp2vtp(self, f, **kwargs)
        return

    def to_shapefile(self, fstem):
        """ Save line to a shapefile """
        if self.rank == 2:
            shp.write_multipoint2(self, fstem)
        elif self.rank == 3:
            shp.write_multipoint3(self, fstem)
        else:
            raise IOError("Rank must be 2 or 3 to write as a shapefile")
        return


class Multipoint(MultipointBase):
    """ Point cloud with associated attributes. This is a base class for the
    polyline and polygon classes.

    *coords*        List of 2-tuples or 3-tuples
    *data*          List, dictionary, or Metadata object for point-specific
                    data [default None]
    *properties*    Dictionary of geometry specific data [default None]
    *crs*           Coordinate reference system instance [default CARTESIAN]
    """
    _geotype = "Multipoint"

    @property
    def __geo_interface__(self):
        return {"type" : "MultiPoint", "bbox" : self.bbox, "coordinates" : self.vertices}

    def within_radius(self, pt, radius):
        """ Return Multipoint of subset that is within *radius* of *pt*.
        """
        distances = self.distances_to(pt)
        indices = [i for i,d in enumerate(distances) if d <= radius]
        return self._subset(indices)

    def within_bbox(self, bbox):
        """ Return Multipoint subset that is within a square bounding box
        given by (xmin, xymin, xmax, ymax). """
        filtbbox = lambda pt: (bbox[0] <= pt.vertex[0] <= bbox[2]) and \
                              (bbox[1] <= pt.vertex[1] <= bbox[3])
        indices = [i for (i, pt) in enumerate(self) if filtbbox(pt)]
        return self._subset(indices)


class ConnectedMultipoint(MultipointBase):
    """ Class for Multipoints in which vertices are assumed to be connected. """

    @property
    def length(self):
        """ Returns the length of the line/boundary. """
        points = [Point(v, crs=self._crs) for v in self.vertices]
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
        interxbool = (np.nan in _vecgeo.intersection(a[0][0], a[1][0], b[0][0], b[1][0],
                                                     a[0][1], a[1][1], b[0][1], b[1][1])
                    for a in self.segments for b in other.segments)
        if self._bbox_overlap(other) and (True not in interxbool):
            return True
        else:
            return False

    def intersections(self, other, keep_duplicates=False):
        """ Return the intersections with another geometry as a Multipoint. """
        interx = (_vecgeo.intersection(a[0][0], a[1][0], b[0][0], b[1][0],
                                          a[0][1], a[1][1], b[0][1], b[1][1])
                     for a in self.segments for b in other.segments)
        if not keep_duplicates:
            interx = set(interx)
        interx_points = []
        for vertex in interx:
            if np.nan not in vertex:
                interx_points.append(Point(vertex, properties=self.properties,
                                           crs=self._crs))
        nandata = [np.nan for _ in interx_points]
        keys = [k for k in set(list(self.data.keys()) + list(other.data.keys()))]
        d = Metadata(dict((key,nandata) for key in keys))
        return Multipoint(interx_points, data=d)

    def shortest_distance_to(self, pt):
        """ Return the shortest distance from any position on the Multipoint
        boundary to *pt* (Point).
        """
        if self._crs != pt._crs:
            raise CRSError("Point must have matching CRS")
        elif self._crs == Cartesian:
            func = _vecgeo.pt_nearest_planar
        else:
            fwd = self._crs.forward
            inv = self._crs.inverse
            func = lambda *args: _vecgeo.pt_nearest_proj(fwd, inv, *args, tol=0.01)

        segments = list(self.segments)
        point_dist = map(func,
                         (pt.vertex for _ in segments),
                         (seg[0].vertex for seg in segments),
                         (seg[1].vertex for seg in segments))
        distances = [i[1] for i in point_dist]
        return min(distances)

    def nearest_on_boundary(self, pt):
        """ Returns the position on the Multipoint boundary that is nearest to
        pt (Point). If two points are equidistant, only one will be returned.
        """
        if self._crs != pt._crs:
            raise CRSError("Point must have matching CRS")
        elif self._crs == Cartesian:
            func = _vecgeo.pt_nearest_planar
        else:
            fwd = self._crs.forward
            inv = self._crs.inverse
            func = lambda *args: _vecgeo.pt_nearest_proj(fwd, inv, *args, tol=0.01)

        segments = list(self.segments)
        point_dist = list(map(func,
                              (pt.vertex for _ in segments),
                              (seg[0].vertex for seg in segments),
                              (seg[1].vertex for seg in segments)))
        distances = [i[1] for i in point_dist]
        imin = distances.index(min(distances))
        return Point(point_dist[imin][0], crs=self._crs)

    def within_distance(self, pt, distance):
        """ Test whether a point is within *distance* of a ConnectedMultipoint. """
        return all(distance >= seg.shortest_distance_to(pt) for seg in self.segments)


class Line(ConnectedMultipoint):
    """ Line composed of connected vertices.

    *coords*        List of 2-tuples or 3-tuples
    *data*          List, dictionary, or Metadata object for point-specific
                    data [default None]
    *properties*    Dictionary of geometry specific data [default None]
    *crs*           Coordinate reference system instance [default CARTESIAN]
    """
    _geotype = "Line"

    @property
    def __geo_interface__(self):
        return {"type" : "LineString", "bbox" : self.bbox, "coordinates" : self.vertices}

    def add_vertex(self, vertex):
        """ Add a vertex to self.vertices. """
        if isinstance(vertex, Point):
            if self.rank == 2:
                self.vertices.append((vertex.x, vertex.y))
            elif self.rank == 3:
                self.vertices.append((vertex.x, vertex.y, vertex.z))
            self.data += vertex.data
        else:
            if self.rank == 2:
                self.vertices.append((vertex[0], vertex[1]))
            elif self.rank == 3:
                self.vertices.append((vertex[0], vertex[1], vertex[2]))
        return self

    def remove_vertex(self, index):
        """ Removes a vertex from the register by index. """
        pt = Point(self.vertices.pop(index), data=self.data[index])
        del self.data[index]
        return pt

    def extend(self, other):
        """ Combine two lines, provided that that the data formats are similar.
        """
        if self.rank == other.rank:
            if self._geotype == other._geotype:
                self.vertices.extend(other.vertices)
                self.data = self.data + other.data
            else:
                GGeoError('Cannot add inconsistent geometry types')
        else:
            GGeoError('Cannot add geometries with inconsistent rank')
        return self

    def cumlength(self):
        """ Returns the cumulative length by segment, prefixed by zero. """
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
        Ltotal = self.cumlength()[-1]
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

    #def direction(self):
    #    """ Returns a vector of the azimuth along a line at each point,
    #    weighted by neighbour position. """
    #    # The weights are calculated to give greater weight to azimuths to
    #    # nearby points
    #    # a' = w1 a1 + w2 a2
    #    #    = (a1 h2) / (h1 + h2) + (a2 h1) / (h1 + h2)
    #    #    = (a1 h2 + a2 h1) / (h1 + h2)
    #    azs = np.empty(len(self)-1)
    #    hs = np.empty(len(self)-1)
    #    pts = [pt for pt in self]
    #    for i in range(len(self)-1):
    #        azs[i] = pt[i].azimuth(pt[i+1])
    #        hs[i] = pt[i].distance(pt[i+1])
    #    weighted_azs = np.empty(len(self))
    #    weighted_azs[0] = azs[0]
    #    weighted_azs[-1] = azs[-1]
    #    for i in range(1, len(self)-1):
    #        raise NotImplementedError("need to phase unwrap")
    #        weighted_azs[i] = (azs[i]*hs[i+1] + azs[i+1]*hs[i]) / (hs[i] + hs[i+1])
    #    return weighted_azs

    def to_polygon(self):
        """ Returns a polygon. """
        return Polygon(self.vertices, data=self.data, properties=self.properties,
                       crs=self._crs)

    def to_shapefile(self, fstem):
        """ Save line to a shapefile """
        if self.rank == 2:
            shp.write_line2(self, fstem)
        elif self.rank == 3:
            shp.write_line3(self, fstem)
        else:
            raise IOError("rank must be 2 or 3 to write as a shapefile")
        return


class Polygon(ConnectedMultipoint):
    """ Polygon, composed of a closed sequence of vertices.

    *coords*        List of 2-tuples or 3-tuples
    *data*          List, dictionary, or Metadata object for point-specific
                    data [default None]
    *properties*    Dictionary of geometry specific data [default None]
    *subs*          List of sub-polygons [default None]
    *crs*           Coordinate reference system instance [default CARTESIAN]
    """
    _geotype = "Polygon"
    subs = []

    def __init__(self, vertices, data=None, properties=None, subs=None, **kwargs):
        vertices = list(vertices)
        ConnectedMultipoint.__init__(self, vertices, data=data,
                                     properties=properties, **kwargs)
        self.subs = subs if subs is not None else []
        return

    def __getitem__(self, key):
        if isinstance(key, slice):
            ind = key.indices(len(self))
            if len(self) != ((ind[1] - ind[0]) // ind[2]):
                return Line(self.vertices[key], data=self.data[key],
                            properties=self.properties, crs=self._crs)
        return super(Polygon, self).__getitem__(key)

    @property
    def __geo_interface__(self):
        coords = [self.vertices]
        for geom in self.subs:
            coords.append(geom.vertices)
        return {"type" : "Polygon", "bbox" : self.bbox, "coordinates" : coords}

    def _subset(self, idxs):
        """ Return a subset defined by index in *idxs*. """
        vertices = [self.vertices[i] for i in idxs]
        data = self.data.sub(idxs)
        subset = Line(vertices, data=data, properties=self.properties,
                      crs=self._crs, copy_metadata=False)
        return subset

    def isclockwise(self):
        """ Return whether polygon winds clockwise around its interior. """
        s = sum((seg[1][0] - seg[0][0]) * (seg[1][1] + seg[0][1])
                for seg in self.segments)
        return s > 0

    @property
    def segments(self):
        """ Returns an generator of adjacent line segments.
        Unique to Polygon: appends a final segment to close the Polygon.
        """
        L = len(self.vertices)
        return itertools.chain((self._subset((i,i+1)) for i in range(len(self)-1)),
                               (self._subset((L-1,0)),))

    @property
    def length(self):
        return self.perimeter

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
        x, y = self.coordinates
        x0 = np.min(x)
        a = (0.5*(x[0] + x[-1]) - x0) * (y[0] - y[-1])
        a += sum((0.5*(x[i+1]+x[i]) - x0) * (y[i+1] - y[i]) for i in range(len(x)-1))
        return abs(a) - sum(map(lambda p: p.area, self.subs))

    @staticmethod
    def _signcross(a, b):
        """ Return sign of 2D cross product a x b """
        c = (a[0]*b[1]) - (a[1]*b[0])
        if c != 0:
            return c/abs(c)
        else:
            return 0

    def _contains_cartesian(self, pt):
        """ Fast `contains` for 2D data on cartesian plane """
        a = np.array(self.vertices[-1])
        b = np.array(self.vertices[0])
        p = np.array(pt.vertex)
        s = self._signcross(b-p, b-a)
        i = 0
        while True:
            a = np.array(self.vertices[i])
            b = np.array(self.vertices[i+1])
            s_ = self._signcross(b-p, b-a)
            i += 1
            if s_ == -s:
                # Point is outside
                return False
            elif i == len(self.vertices)-1:
                # Point in inside, so check subpolys
                return not any(self._contains_cartesiah(p, pt) for p in self.subs)

    def _contains_projected(self, pt):
        """ Slower `contains` for projected data based on counting
        intersections. """
        nintx = 0
        x, y = pt[0], pt[1]
        for seg in self.segments:
            (a, b) = seg
            if _vecgeo.intersects_cn(x, y, a[0], b[0], a[1], b[1]):
                nintx += 1

        # If odd, point is inside so check subpolys
        return nintx % 2 == 1 and not any(p.contains(pt) for p in self.subs)

    def contains(self, pt):
        """ Returns True if pt is inside or on the boundary of the
        polygon, and False otherwise. Uses a crossing number scheme.
        """
        if self.crs is Cartesian:
            return self._contains_cartesian(pt)
        else:
            return self._contains_projected(pt)

    def to_polyline(self):
        """ Returns a self-closing polyline. Discards sub-polygons. """
        return Line(self.vertices, properties=self.properties, data=self.data,
                    crs=self.crs)

    def to_shapefile(self, fstem):
        """ Save line to a shapefile """
        if self.rank == 2:
            shp.write_poly2(self, fstem)
        elif self.rank == 3:
            shp.write_poly3(self, fstem)
        else:
            raise IOError("rank must be 2 or 3 to write as a shapefile")
        return


class GeometryError(Exception):
    """ Base class for geometry module errors. """
    def __init__(self, message=''):
        self.message = message
    def __str__(self):
        return self.message


class GInitError(GeometryError):
    """ Exception to raise when a geometry object fails to initialize. """
    def __init__(self, message=''):
        self.message = message


class GUnitError(GeometryError):
    """ Exception to raise there is a projected unit problem. """
    def __init__(self, message=''):
        self.message = message


class GGeoError(GeometryError):
    """ Exception to raise when a geometry object attempts an invalid transform. """
    def __init__(self, message=''):
        self.message = message


def points_to_multipoint(points):
    """ Merge *points* into a Multipoint instance. Point properties are stored
    as Multipoint data. All points must use the same CRS.
    """
    crs = points[0]._crs
    if not all(pt._crs == crs for pt in points):
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

