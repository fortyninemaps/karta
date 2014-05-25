"""
Geographical measurement and simple analysis module for Python 2.X.X. Provides
Point, Multipoint, Line, and Polygon classes, with methods for simple
measurements such as distance, area, and bearing.

Written by Nat Wilson (njwilson23@gmail.com)
"""

import math
import sys
import traceback
import numpy as np
import shapefile
from . import vtk
from . import geojson
from . import xyfile
from . import _shpfuncs
from .. import crs as kcrs

from collections import deque
from .metadata import Metadata

try:
    from . import _cvectorgeo as _vecgeo
except ImportError:
    sys.stderr.write("falling back on slow _vectorgeo")
    from . import _vectorgeo as _vecgeo

try:
    import shapely.geometry as geometry
except ImportError:
    pass

try:
    import pyproj
    PYPROJ = True
    geod = pyproj.Geod(ellps="WGS84")
except ImportError:
    PYPROJ = False

class Geometry(object):
    """ This is the abstract base class for all geometry types, i.e. Point,
    Multipoints and subclasses thereof. """
    _geotype = None

    def __init__(self, crs=kcrs.CARTESIAN):
        self.properties = {}
        self._crs = crs
        return

    def _distance(self, pos0, pos1):
        """ Generic method for calculating distance between positions that
        respects CRS """
        if self._crs == kcrs.CARTESIAN:
            dist = _vecgeo.distance((pos0.x, pos0.y), (pos1.x, pos1.y))
        elif PYPROJ:
            _, _, dist = geod.inv(pos0.x, pos0.y, pos1.x, pos1.y, radians=False)
        else:
            dist = greatcircle(pos0, pos1)
        return dist

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
    """ This defines the point class, from which x,y[,z] points can be
    constructed. """
    _geotype = "Point"

    def __init__(self, coords, data=None, properties=None, **kwargs):
        if properties is None: properties = {}
        if not hasattr(coords, "__iter__"):
            raise ValueError("input argument must be a sequence")
        super(Point, self).__init__(**kwargs)
        self.vertex = coords
        self._setxyz()
        self.data = Metadata(data, singleton=True)
        if hasattr(properties, "keys"):
            self.properties = properties
        return

    def __getitem__(self, idx):
        return self.vertex[idx]

    def __repr__(self):
        return 'Point(' + str(self.vertex) + ')'

    def __eq__(self, other):
        if hasattr(other, "vertex"):
            return (tuple(self.vertex) == tuple(other.vertex)) and \
                   (self.data == other.data) and \
                   (self.properties == other.properties)
        else:
            return False

    def _setxyz(self):
        """ Create convenience *x*, *y* (, *z*) attributes from vertices. """
        self.x = self.vertex[0]
        self.y = self.vertex[1]
        try:
            self.z = self.vertex[2]
            self.rank = 3
        except IndexError:
            self.z = None
            self.rank = 2
        return

    def isat(self, other, tol=0.001):
        """ Test whether point vertices are the same, regardless of properties or data.
        Tolarance is not geodetic. """
        return False not in [abs(a-b) <= tol for (a, b) in zip(self.vertices,
                                                               other.vertices)]

    def get_vertex(self):
        """ Return the Point vertex as a tuple. """
        return self.vertex

    def coordsxy(self, convert_to=False):
        """ Returns the x,y coordinates. Convert_to may be set to 'deg'
        or 'rad' for convenience.  """
        if convert_to == 'rad':
            return (self.x*3.14159/180., self.y*3.14159/180.)
        elif convert_to == 'deg':
            return (self.x/3.14159*180., self.y/3.14159*180.)
        else:
            return (self.x, self.y)

    def azimuth(self, other, spherical=False):
        """ Returns the compass azimuth from self to other in radians (i.e.
        clockwise, with north at 0). Returns NaN if points are coincident. """

        if self.coordsxy() == other.coordsxy():
            return np.nan

        if self._crs == kcrs.LONLAT and other._crs == kcrs.LONLAT:
            if not PYPROJ:
                raise CRSError("Non-cartesian points require pyproj")
            az1, _, _ = geod.inv(self.x, self.y, other.x, other.y)
            return az1 * math.pi / 180.0

        elif self._crs == kcrs.CARTESIAN and other._crs == kcrs.CARTESIAN:
            dx = other.x - self.x
            dy = other.y - self.y

            if dx > 0:
                if dy > 0:
                    return math.atan(dx/dy)
                elif dy < 0:
                    return math.pi - math.atan(-dx/dy)
                else:
                    return 0.5*math.pi
            elif dx < 0:
                if dy > 0:
                    return 2*math.pi - math.atan(-dx/dy)
                elif dy < 0:
                    return math.pi + math.atan(dx/dy)
                else:
                    return 1.5*math.pi
            else:
                if dy > 0:
                    return 0.0
                else:
                    return math.pi

        else:
            raise CRSError("Azimuth undefined for points in CRS {0} and "
                           "{1}".format(self._crs, other._crs))

    def walk(self, distance, direction):
        """ Returns the point reached when moving in a given direction for
        a given distance from a specified starting location.

            distance (float): distance to walk
            direction (float): horizontal walk direction in radians
        """
        if self._crs == kcrs.CARTESIAN:
            dx = distance * math.cos(direction)
            dy = distance * math.sin(direction)

            if self.rank == 2:
                return Point((self.x+dx, self.y+dy),
                             properties=self.properties, data=self.data, crs=self._crs)
            elif self.rank == 3:
                return Point((self.x+dx, self.y+dy, self.z),
                             properties=self.properties, data=self.data, crs=self._crs)

        elif PYPROJ:
            (x, y, backaz) = geod.fwd(self.x, self.y, direction*180.0/math.pi,
                                      distance, radians=False)
            return Point((x, y), properties=self.properties, data=self.data,
                         crs=self._crs)

        else:
            raise CRSError("Non-cartesian points require pyproj")

    def distance(self, other):
        """ Returns a distance to another Point. If the coordinate system is
        geographical and a third (z) coordinate exists, it is assumed to have
        the same units as the real-world horizontal distance (i.e. meters). """
        if self._crs != other._crs:
            raise CRSError("Points must share the same coordinate system.")
        flat_dist = self._distance(self, other)
        if None in (self.z, other.z):
            return flat_dist
        else:
            return math.sqrt(flat_dist**2. + (self.z-other.z)**2.)

    def greatcircle(self, other):
        """ Return the great circle distance between two geographical points. """
        if not PYPROJ:
            raise CRSError("Non-cartesian points require pyproj")
        if not (self._crs == kcrs.LONLAT and other._crs == kcrs.LONLAT):
            raise CRSError("Great circle distances require both points to be "
                           "in geographical coordinates")
        az1, az2, dist = geod.inv(self.x, self.y, other.x, other.y, radians=False)
        return dist

    def shift(self, shift_vector):
        """ Shift point by the amount given by a vector. Operation occurs
        in-place """
        if len(shift_vector) != self.rank:
            raise GGeoError('Shift vector length must equal geometry rank.')

        self.vertex = tuple([a+b for a,b in zip(self.vertex, shift_vector)])
        self._setxyz()
        return self

    def as_geojson(self, **kwargs):
        """ Write data as a GeoJSON string to a file-like object `f`.

        Parameters
        ----------
        f : file-like object to recieve the GeoJSON string

        *kwargs* include:
        crs : coordinate reference system
        crs_fmt : format of `crs`; may be one of ('epsg','ogc_crs_urn')
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
        crs_fmt : format of `crs`; may be one of ('epsg','ogc_crs_urn')
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

    def to_shapely(self):
        """ Returns a Shapely Point instance. """
        try:
            return geometry.Point(self.x, self.y, self.z)
        except NameError:
            raise ImportError('Shapely module did not import\n')


class MultipointBase(Geometry):
    """ Point cloud with associated attributes. This is a base class for the
    polyline and polygon classes. """
    _geotype = "MultipointBase"

    def __init__(self, vertices, data=None, properties=None, **kwargs):
        """ Create a feature with multiple vertices.

        vertices : a list of tuples containing point coordinates.

        data : is either `None` a list of point attributes, or a dictionary of
        point attributes. If `data` is not `None`, then it (or its values) must
        match `vertices` in length.
        """
        super(MultipointBase, self).__init__(**kwargs)
        vertices = list(vertices)
        if properties is None:
            properties = {}
        if len(vertices) > 0:
            self.rank = len(vertices[0])

            if self.rank > 3 or self.rank < 2:
                raise GInitError('Input must be doubles or triples\n')
            elif False in [self.rank == len(i) for i in vertices]:
                raise GInitError('Input must have consistent rank\n')
            else:
                self.vertices = [tuple(i) for i in vertices]

        else:
            self.rank = None
            self.vertices = []

        if hasattr(properties, "keys"):
            self.properties = properties

        self.data = Metadata(data)
        return

    def __repr__(self):
        if len(self) < 5:
            ppverts = str(self.vertices)
        else:
            ppverts = str(self.vertices[:2])[:-1] + "..." + str(self.vertices[-2:])[1:]
        return '{typ}({verts})>'.format(
                typ=str(type(self))[:-1], verts=ppverts, prop=self.properties)

    def __str__(self):
        if len(self) < 5:
            ppverts = str(self.vertices)
        else:
            ppverts = str(self.vertices[:2])[:-1] + "..." + str(self.vertices[-2:])[1:]
        return '{typ}\n{verts}\nProperties:{prop}'.format(
                typ=type(self), verts=ppverts, prop=self.properties)

    def __len__(self):
        return len(self.vertices)

    def __getitem__(self, key):
        if isinstance(key, (int, np.int64)):
            return Point(self.vertices[key], data=self.data[key],
                         properties=self.properties, crs=self._crs)
        elif isinstance(key, slice):
            return type(self)(self.vertices[key], data=self.data[key],
                              properties=self.properties, crs=self._crs)
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
               (False not in (a==b for a,b in zip(self, other)))

    def _bbox_overlap(self, other):
        """ Return whether bounding boxes between self and another geometry
        overlap.
        """
        reg0 = self.get_bbox()
        reg1 = other.get_bbox()
        return (reg0[0] < reg1[1] and reg0[1] > reg1[0] and
                reg0[2] < reg1[3] and reg0[3] > reg1[2])

    def get_bbox(self):
        """ Return the extents of a bounding box as
            (xmin, ymax, ymin, ymax, [zmin, zmin]).
        """
        if self.rank == 2:
            x, y = self.get_coordinate_lists()
            bbox = (min(x), max(x), min(y), max(y))
        elif self.rank == 3:
            x, y, z = self.get_coordinate_lists()
            bbox = (min(x), max(x), min(y), max(y), min(z), max(z))
        return bbox

    def print_vertices(self):
        """ Prints an enumerated list of indices. """
        for i, vertex in enumerate(self.vertices):
            print("{0}\t{1}".format(i, vertex))

    def get_vertices(self):
        """ Return vertices as a list of tuples. """
        return np.array(self.vertices)

    #def get_data(self, fields=None):
    #    """ Return data as an array, regardless of internal type. Optionally
    #    takes the keyword argument *fields*, which is an iterable listing the
    #    columns from the data dictionary to retrieve. """
    #    if hasattr(self.data, 'keys') and hasattr(self.data.values, '__call__'):
    #        if fields is not None:
    #            data = np.array([self.data[key] for key in fields])
    #        else:
    #            data = np.array(self.data.values())
    #    else:
    #        data = np.array(self.data)
    #    return data.T

    def get_coordinate_lists(self):
        """ Return X, Y, and Z lists. If self.rank == 2, Z will be
        zero-filled. """
        X = [i[0] for i in self.vertices]
        Y = [i[1] for i in self.vertices]
        if self.rank == 3:
            Z = [i[2] for i in self.vertices]
            return X, Y, Z
        else:
            return X, Y

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

    def _subset(self, idxs):
        """ Return a subset defined by index in *idxs*. """
        vertices = [self.vertices[i] for i in idxs]
        data = self.data.sub(idxs)
        subset = type(self)(vertices, data=data, properties=self.properties, crs=self._crs)
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

    def greatcircles_to(self, pt):
        """ Return the great circle distances from each vertex to a point. """
        d = [pt.greatcircle(a) for a in self]
        return np.array(d)

    def nearest_point_to(self, pt):
        """ Returns the internal point that is nearest to pt (Point class).
        If two points are equidistant, only one will be returned.
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

    # This code should compute the convex hull of the points and then test the
    # hull's combination space
    # Alternatively, calculate the eigenvectors, rotate, and cherrypick the
    # points
    #def max_dimension(self):
    #    """ Return the two points in the Multipoint that are furthest
    #    from each other. """
    #    dist = lambda xy0, xy1: math.sqrt((xy1[0]-xy0[0])**2 +
    #                                      (xy1[1]-xy0[1])**2)

    #    P = [(p0, p1) for p0 in self.vertices for p1 in self.vertices]
    #    D = map(dist, (p[0] for p in P), (p[1] for p in P))
    #    return P[D.index(max(D))]

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
        crs_fmt : format of `crs`; may be one of ('epsg','ogc_crs_urn')
        bbox : an optional bounding box tuple in the form (w,e,s,n)
        """
        writer = geojson.GeoJSONWriter(self, crs=self._crs,
                    crs_fmt="ogc_crs_urn", **kwargs)
        return writer.print_json()

    def to_geojson(self, f, **kwargs):
        """ Write data as a GeoJSON string to a file-like object `f`.

        Parameters
        ----------
        f : file-like object to recieve the GeoJSON string

        *kwargs* include:
        crs : coordinate reference system
        crs_fmt : format of `crs`; may be one of ('epsg','ogc_crs_urn')
        bbox : an optional bounding box tuple in the form (w,e,s,n)
        """
        try:
            if not hasattr(f, "write"):
                fobj = open(f, "w")
            else:
                fobj = f
            writer = geojson.GeoJSONWriter(self, crs=self._crs,
                        crs_fmt="ogc_crs_urn", **kwargs)
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
            _shpfuncs.write_multipoint2(self, fstem)
        elif self.rank == 3:
            _shpfuncs.write_multipoint3(self, fstem)
        else:
            raise IOError("rank must be 2 or 3 to write as a shapefile\n")
        return

class Multipoint(MultipointBase):
    """ Point cloud with associated attributes. This is a base class for the
    polyline and polygon classes. """
    _geotype = "Multipoint"

    def __init__(self, vertices, data=None, properties=None, **kwargs):
        """ Create a feature with multiple vertices.

        vertices : a list of tuples containing point coordinates.

        data : is either `None` a list of point attributes, or a dictionary of
        point attributes. If `data` is not `None`, then it (or its values) must
        match `vertices` in length.
        """
        if False not in (hasattr(v, "_geotype") and \
                         v._geotype == "Point" and \
                         v.rank == vertices[0].rank for v in vertices):
            points = vertices
            crs = points[0]._crs
            if False in (pt._crs == crs for pt in points):
                raise CRSError("All points must share the same CRS")

            keys = list(points[0].properties.keys())
            for pt in points[1:]:
                for key in keys:
                    if key not in pt.properties:
                        keys.pop(keys.index(key))

            ptdata = {}
            for key in keys:
                ptdata[key] = [pt.properties[key] for pt in points]

            if data is not None:
                ptdata.update(data)
            self.data = Metadata(ptdata)
            self.vertices = [pt.vertex for pt in points]
            self.rank = vertices[0].rank
            self.properties = properties
            self._crs = crs

        else:

            super(Multipoint, self).__init__(vertices,
                                             data=data, 
                                             properties=properties,
                                             **kwargs)
        return

    def within_radius(self, pt, radius):
        """ Return Multipoint of subset that is within *radius* of *pt*.
        """
        distances = self.distances_to(pt)
        indices = [i for i,d in enumerate(distances) if d <= radius]
        return self._subset(indices)

    def within_bbox(self, bbox):
        """ Return Multipoint subset that is within a square boundaing box
        given by (xmin, xmax, ymin, ymax). """
        filtbbox = lambda pt: (bbox[0] <= pt.vertex[0] <= bbox[1]) and \
                              (bbox[2] <= pt.vertex[1] <= bbox[3])
        indices = [i for (i, pt) in enumerate(self) if filtbbox(pt)]
        return self._subset(indices)


class ConnectedMultipoint(MultipointBase):
    """ Class for Multipoints in which vertices are assumed to be connected. """

    def length(self):
        """ Returns the length of the line/boundary. """
        points = [Point(i, crs=self._crs) for i in self.vertices]
        distances = [a.distance(b) for a, b in zip(points[:-1], points[1:])]
        return sum(distances)

    def segments(self):
        """ Returns an iterator of adjacent line segments. """
        return (self._subset((i,i+1)) for i in range(len(self)-1))

    def segment_tuples(self):
        """ Returns an iterator of adjacent line segments as coordinate tuples. """
        return ((self.vertices[i], self.vertices[i+1])
                for i in range(len(self.vertices)-1))

    def intersects(self, other):
        """ Return whether an intersection exists with another geometry. """
        interxbool = (_vecgeo.intersects(a[0][0], a[1][0], b[0][0], b[1][0],
                                         a[0][1], a[1][1], b[0][1], b[1][1])
                    for a in self.segments() for b in other.segments())
        if self._bbox_overlap(other) and (True in interxbool):
            return True
        else:
            return False

    def intersections(self, other):
        """ Return the intersections with another geometry. """
        interx = (_vecgeo.intersections(a[0][0], a[1][0], b[0][0], b[1][0],
                                        a[0][1], a[1][1], b[0][1], b[1][1])
                    for a in self.segments() for b in other.segments())
        return list(filter(lambda a: np.nan not in a, interx))

    def shortest_distance_to(self, pt):
        """ Return the shortest distance from a point on the Multipoint
        boundary to *pt* (Point) """
        point_dist = map(_vecgeo.pt_nearest,
                                [pt.vertex for seg in self.segments()],
                                [seg[0].vertex for seg in self.segments()],
                                [seg[1].vertex for seg in self.segments()])
        distances = [i[1] for i in point_dist]
        return min(distances)

    def nearest_on_boundary(self, pt):
        """ Returns the point on the Multipoint boundary that is nearest to pt
        (Point).

        Warning: If two points are equidistant, only one will be returned.
        """
        point_dist = list(map(_vecgeo.pt_nearest,
                                [pt.vertex for seg in self.segments()],
                                [seg[0].vertex for seg in self.segments()],
                                [seg[1].vertex for seg in self.segments()]))
        distances = [i[1] for i in point_dist]
        return Point(point_dist[distances.index(min(distances))][0],
                     properties=self.properties, crs=self._crs)

    def within_distance(self, pt, distance):
        """ Test whether a point is within *distance* of a ConnectedMultipoint. """
        #return distance <= self.shortest_distance_to(pt)
        return True in (distance >= seg.shortest_distance_to(pt)
                        for seg in self.segments())

class Line(ConnectedMultipoint):
    """ This defines the polyline class, from which geographic line
    objects can be constructed. Line objects consist of joined,
    georeferenced line segments.
    """
    _geotype = "Line"

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

    def displacement(self):
        """ Returns the distance between the first and last vertex. """
        return Point(self.vertices[0]).distance(Point(self.vertices[-1]))

    def direction(self):
        """ Returns a vector of the azimuth along a line at each point,
        weighted by neighbour position. """
        # The weights are calculated to give greater weight to azimuths to
        # nearby points
        # a' = w1 a1 + w2 a2
        #    = (a1 h2) / (h1 + h2) + (a2 h1) / (h1 + h2)
        #    = (a1 h2 + a2 h1) / (h1 + h2)
        azs = np.empty(len(self)-1)
        hs = np.empty(len(self)-1)
        pts = [pt for pt in self]
        for i in range(len(self)-1):
            azs[i] = pt[i].azimuth(pt[i+1])
            hs[i] = pt[i].distance(pt[i+1])
        weighted_azs = np.empty(len(self))
        weighted_azs[0] = azs[0]
        weighted_azs[-1] = azs[-1]
        for i in range(1, len(self)-1):
            raise NotImplementedError("need to phase unwrap")
            weighted_azs[i] = (azs[i]*hs[i+1] + azs[i+1]*hs[i]) / (hs[i] + hs[i+1])
        return weighted_azs

    def to_polygon(self):
        """ Returns a polygon. """
        return Polygon(self.vertices)

    def to_shapely(self):
        """ Returns a Shapely LineString instance. """
        try:
            if self.rank == 2:
                return geometry.LineString([(v[0], v[1]) for v in self.vertices])
            elif self.rank == 3:
                return geometry.LineString([(v[0], v[1], v[2]) for v in self.vertices])
        except NameError:
            raise GuppyError('Shapely module not available\n')

    def to_shapefile(self, fstem):
        """ Save line to a shapefile """
        if self.rank == 2:
            _shpfuncs.write_line2(self, fstem)
        elif self.rank == 3:
            _shpfuncs.write_line3(self, fstem)
        else:
            raise IOError("rank must be 2 or 3 to write as a shapefile\n")
        return

class Polygon(ConnectedMultipoint):
    """ This defines the polygon class, from which geographic
    polygons objects can be created. Polygon objects consist of
    point nodes enclosing an area.
    """
    _geotype = "Polygon"
    subs = []

    def __init__(self, vertices, data=None, properties=None, subs=None, **kwargs):
        ConnectedMultipoint.__init__(self, vertices, data=data, properties=properties, **kwargs)
        if vertices[0] != vertices[-1]:
            self.vertices.append(vertices[0])
        self.subs = subs if subs is not None else []
        return

    def __getitem__(self, key):
        if isinstance(key, slice):
            ind = key.indices(len(self))
            if len(self) != ((ind[1] - ind[0]) // ind[2]):
                return Line(self.vertices[key], data=self.data[key], 
                            properties=self.properties, crs=self._crs)
        return super(Polygon, self).__getitem__(key)

    def perimeter(self):
        """ Return the perimeter of the polygon. If there are sub-polygons,
        their perimeters are added recursively. """
        return self.length() + sum([p.perimeter() for p in self.subs])

    def area(self):
        """ Return the two-dimensional area of the polygon. If there are
        sub-polygons, their areas are subtracted. """
        a = 0.0
        for i in range(len(self.vertices)-1):
            a += 0.5 * abs((self.vertices[i][0] + self.vertices[i+1][0])
                         * (self.vertices[i][1] - self.vertices[i+1][1]))
        return a - sum(map(lambda p: p.area(), self.subs))

    def contains(self, pt):
        """ Returns True if pt is inside or on the boundary of the
        polygon, and False otherwise.
        """
        def possible(pt, v1, v2):
            """ Quickly assess potential for an intersection with an x+
            pointing ray. """
            x = pt.vertex[0]
            y = pt.vertex[1]
            if ( ((y > v1[1]) is not (y > v2[1]))
            and ((x < v1[0]) or (x < v2[0])) ):
                return True
            else:
                return False

        # Find how many boundaries a ray pointing out from point crosses
        bool2int = lambda tf: (tf and 1 or 0)
        rvertices = deque(self.vertices)
        rvertices.rotate(1)
        segments = [(v1, v2) for v1, v2 in zip(self.vertices, rvertices)
                    if possible(pt, v1, v2)]

        n_intersect = sum([bool2int(
            isinstance(ray_intersection(pt.vertex, seg[0], seg[1]), tuple))
            for seg in segments])

        if n_intersect % 2 == 1:    # If odd, point is inside so check subpolys
            if True not in (p.contains(pt) for p in self.subs):
                return True
        return False                # Point was outside or was in a subpoly

    def to_polyline(self):
        """ Returns a self-closing polyline. Discards sub-polygons. """
        return Line(self.vertices)

    def to_shapely(self):
        """ Returns a Shapely Polygon instance. """
        try:
            shp = geometry.Polygon(self.vertices,
                                   interiors=[p.vertices for p in self.subs])
        except NameError:
            raise ImportError('Shapely module did not import\n')
        return shp

    def to_shapefile(self, fstem):
        """ Save line to a shapefile """
        if self.rank == 2:
            _shpfuncs.write_poly2(self, fstem)
        elif self.rank == 3:
            _shpfuncs.write_poly3(self, fstem)
        else:
            raise IOError("rank must be 2 or 3 to write as a shapefile\n")
        return


class GuppyError(Exception):
    """ Base class for guppy module errors. """
    def __init__(self, message=''):
        self.message = message
    def __str__(self):
        return self.message


class GInitError(GuppyError):
    """ Exception to raise when a guppy object fails to initialize. """
    def __init__(self, message=''):
        self.message = message


class GUnitError(GuppyError):
    """ Exception to raise there is a projected unit problem. """
    def __init__(self, message=''):
        self.message = message


class GGeoError(GuppyError):
    """ Exception to raise when a guppy object attempts an invalid transform. """
    def __init__(self, message=''):
        self.message = message

class CRSError(GuppyError):
    """ Exception to raise for invalid geodetic operations. """
    def __init__(self, message=''):
        self.message = message

def ray_intersection(pt, endpt1, endpt2, direction=0.0):
    """ Determines whether a ray intersects a line segment. If yes,
    returns the point of intersection. If no, return None. Input
    "points" should be tuples or similar. Input direction is in
    radians, and defines the ray path from pt.
    """
    m_ray = math.tan(direction)
    if endpt2[0] != endpt1[0]:
        m_lin = float(endpt2[1] - endpt1[1]) / float(endpt2[0] - endpt1[0])
        if m_ray == m_lin:      # Lines are parallel
            return
        else:
            x_int = ( (m_ray * pt[0] - m_lin * endpt1[0] - pt[1] + endpt1[1])
                / (m_ray - m_lin) )

        # Test that y_int is within segment and ray points toward segment
        if ( (x_int >= endpt1[0]) is not (x_int >= endpt2[0]) and
            (x_int-pt[0] > 0) is (math.cos(direction) > 0) ):
            y_int = ( (m_ray * m_lin * pt[0] - m_ray * m_lin * endpt1[0]
                + m_ray * endpt1[1] - m_lin * pt[1]) / (m_ray - m_lin) )
            return (x_int, y_int)
        else:
            return

    else:       # Line segment is vertical
        if direction % math.pi/2. == 0.0 and direction != 0.0:
            # Lines are parallel
            return
        x_int = float(endpt1[0])
        y_int = (x_int - pt[0]) * m_ray + pt[1]
        # Test that y_int is within segment and ray points toward segment
        if ( (y_int >= endpt1[1]) is not (y_int >= endpt2[1]) and
            (x_int-pt[0] > 0) is (math.cos(direction) > 0) ):
            return (x_int, y_int)
        else:
            return


def greatcircle(pta, ptb, method="vicenty"):
    """ Computes the great circle distance between n point pairs on a
    sphere. Returns a list of length (n-1)

    [pntlist] contains a list of point objects

    [method] may be "vicenty" (default) or "haversine". The Haversine
    method is roughly 20% faster, but may yield rounding errors when
    coordinates are antipodal.
    """

    radius = 6371.
    deg2rad = np.pi / 180.

    x1 = pta.x * deg2rad
    x2 = ptb.x * deg2rad
    y1 = pta.y * deg2rad
    y2 = ptb.y * deg2rad

    dx = x2 - x1
    dy = y2 - y1

    if method == "haversine":
        try:
            distance = 2 * radius * math.asin(math.sqrt((math.sin(dy /
                2.))**2 + math.cos(y1) * math.cos(y2) *
                (math.sin(dx / 2.))**2))
        except GGeoError:
            traceback.print_exc()

    elif method == "vicenty":
        try:
            a = math.sqrt((math.cos(y2) * math.sin(dx))**2 +
                (math.cos(y1) * math.sin(y2) - math.sin(y1) *
                math.cos(y2) * math.cos(dx))**2)
            b = (math.sin(y1) * math.sin(y2) + math.cos(y1) *
                math.cos(y2) * math.cos(dx))
            distance = radius * math.atan2(a, b)
        except ZeroDivisionError:
            raise GGeoError("Zero in denominator")
            return None
        except:
            traceback.print_exc()
    else:
        raise Exception("Distance method unrecognized")
        distance = np.nan

    return distance


def sortby(A, B):
    """ Sort a list A by the values in an ordered list B. """
    if len(A) != len(B):
        raise GGeoError("A and B must be of the same length")
    comb = zip(B,A)
    comb.sort()
    return [i[1] for i in comb]


def points_to_multipoint(points):
    """ Merge *points* into a Multipoint instance. Point properties are stored
    as Multipoint data. All points must use the same CRS.
    """
    crs = points[0]._crs
    if False in (pt._crs == crs for pt in points):
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

