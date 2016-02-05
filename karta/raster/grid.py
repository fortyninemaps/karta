""" Classes for basic grid types """

import copy
import math
import numbers
import warnings
import numpy as np
from . import _gtiff
from .. import errors
from ..crs import Cartesian

try:
    from scipy import interpolate
    HASSCIPY = True
except ImportError:
    HASSCIPY = False

IntegerType = (numbers.Integral, np.int32, np.int64)

class Grid(object):
    """ Grid base class. Intended only for implementing subclasses. The primary
    attributes defined by all Grid-derived classed are _hdr and data.
    """

    @property
    def nodata(self):
        return self._nodata

    def max(self):
        """ Return the maximum non-nan in self.data. """
        return np.nanmax(self.values)

    def min(self):
        """ Return the minimum non-nan in self.data. """
        return np.nanmin(self.values)

    def minmax(self):
        """ Return the minimum and maximum value of data array. """
        return (self.min(), self.max())

    def apply(self, func):
        """ Apply a function *func* to grid values """
        g = self.copy()
        g.values = func(g.values)
        return g

    def copy(self):
        """ Return a copy. """
        return copy.deepcopy(self)


class RegularGrid(Grid):
    """ Regular (structured) grid class. A RegularGrid contains a fixed number
    of rows and columns with a constant spacing and a scalar or vector field
    defined as `values`.

    Positions on a RegularGrid are referenced to their array indices *(i,j)*
    using an affine transform *T* such that

    (a, b, c, d, e, f) = T

    X = a + j * c + i * e
    Y = b + i * d + j * f

    in which *(a,b)* defines the lower-left grid cell corner, *(c,d)* defines
    the resolution in the horizontal and vertical directions, and *(e,f)* can
    be used to define a rotation. In the common case of a "north-up" grid,

    e = f = 0
    """
    def __init__(self, transform, values=None, crs=None, nodata_value=None):
        """ Create a RegularGrid instance with cells referenced according to
        *transform*, which is an iterable or dictionary consisting of
        (xllcenter, yllcenter, dx, dy, xrot, yrot).

        Optional parameters:
        --------------------
        values : grid data (nrows x ncols) [default NaN]
        crs : Grid coordinate reference system [default None]
        """
        if hasattr(transform, "keys"):
            self._transform = tuple([float(transform[f]) for f in
                    ("xllcorner", "yllcorner", "dx", "dy", "xrot", "yrot")])
        elif len(transform) == 6:
            self._transform = tuple(float(a) for a in transform)
        elif len(transform) == 4:
            self._transform = (float(transform[0]), float(transform[1]),
                               float(transform[2]), float(transform[3]), 0, 0)
        else:
            raise errors.GridError("RegularGrid must be initialized with a "
                                   " transform iterable or dictionary")

        if values is not None:
            self.values = values
        else:
            self.values = np.atleast_2d([np.nan])

        if crs is None:
            self.crs = Cartesian
        else:
            self.crs = crs

        if nodata_value is None:
            self._nodata = get_nodata(self.values.dtype.type)
        else:
            self._nodata = nodata_value
        return

    def __add__(self, other):
        if self._equivalent_structure(other):
            return RegularGrid(copy.copy(self.transform),
                               values=self.values+other.values,
                               crs=self.crs, nodata_value=self.nodata)
        else:
            raise errors.NonEquivalentGridError(self, other)

    def __sub__(self, other):
        if self._equivalent_structure(other):
            return RegularGrid(copy.copy(self.transform),
                               values=self.values-other.values,
                               crs=self.crs, nodata_value=self.nodata)
        else:
            raise errors.NonEquivalentGridError(self, other)

    def _equivalent_structure(self, other):
        return (self._transform == other._transform) and \
               (self.values.shape == other.values.shape)

    @property
    def transform(self):
        return self._transform

    @property
    def size(self):
        return self.values.shape

    def center_llref(self):
        """ Return the 'lower-left' reference in terms of a center coordinate.
        """
        T = self._transform
        xllcenter = self._transform[0] + 0.5 * (T[2] + T[4])
        yllcenter = self._transform[1] + 0.5 * (T[3] + T[5])
        return xllcenter, yllcenter

    def corner_llref(self):
        """ Return the 'lower-left' reference in terms of a center coordinate.
        """
        return self._transform[:2]

    def coordmesh(self, *args, **kwargs):
        """ Alias for center_coords() """
        return self.center_coords(*args, **kwargs)

    def center_coords(self):
        t = self._transform
        xcoords = np.empty(self.size)
        ycoords = np.empty(self.size)
        icol = np.arange(self.size[1])
        for i in range(self.size[0]):
            xcoords[i,:] = t[0] + (icol+0.5)*t[2] + (i+0.5)*t[4]
            ycoords[i,:] = (t[1] + (i+0.5)*t[3] + (icol+0.5)*t[5])
        return xcoords, ycoords

    def vertex_coords(self):
        """ Return the coordinates of cell vertices. """
        t = self._transform
        ny, nx = self.size
        xcoords = np.empty((ny+1, nx+1))
        ycoords = np.empty((ny+1, nx+1))
        icol = np.arange(0, self.size[1]+1)
        for i in range(self.size[0]):
            xcoords[i,:] = t[0] + icol*t[2] + i*t[4]
            ycoords[i,:] = (t[1] + i*t[3] + icol*t[5])
        return xcoords, ycoords

    @property
    def bbox(self):
        extent = self.get_extent(reference="edge")
        return extent[0], extent[2], extent[1], extent[3]

    @property
    def data_bbox(self):
        extent = self.get_data_extent(reference="edge")
        return extent[0], extent[2], extent[1], extent[3]

    @property
    def extent(self):
        return self.get_extent()

    def get_bbox(self, crs=None):
        a, b, c, d = self.get_extent(reference="edge", crs=crs)
        return a, c, b, d

    def get_extent(self, reference='center', crs=None):
        """ Return the region characteristics as a tuple (xmin, xmax, ymin,
        ymax). *reference* is a string and may be 'center' or 'edge'. """
        if reference == 'center':
            x0, y0 = self.center_llref()
            n = -1
        elif reference == 'edge':
            x0, y0 = self.corner_llref()
            n = 0
        else:
            raise errors.GridError("`reference` must be 'center' or 'edge'")
        ny, nx = self.values.shape[:2]
        dx, dy = self._transform[2:4]
        sx, sy = self._transform[4:]

        sgn = lambda a: 0 if a==0 else a/abs(a)
        if sgn(dx) == sgn(sx):
            x1 = x0 + dx*(nx+n) + sx*(ny+n)
        else:
            x1 = x0 + dx*(nx+n) - sx*(ny+n)

        if sgn(dy) == sgn(sy):
            y1 = y0 + dy*(ny+n) + sy*(nx+n)
        else:
            y1 = y0 + dy*(ny+n) - sy*(nx+n)

        if crs is None or (crs == self.crs):
            (a, b, c, d) = min(x0, x1), max(x0, x1), min(y0, y1), max(y0, y1)
        else:
            xg0, yg0 = crs.project(*self.crs.project(x0, y0, inverse=True))
            xg1, yg1 = crs.project(*self.crs.project(x0, y1, inverse=True))
            xg2, yg2 = crs.project(*self.crs.project(x1, y0, inverse=True))
            xg3, yg3 = crs.project(*self.crs.project(x1, y1, inverse=True))
            a = min(xg0, xg1, xg2, xg3)
            b = max(xg0, xg1, xg2, xg3)
            c = min(yg0, yg1, yg2, yg3)
            d = max(yg0, yg1, yg2, yg3)
        return a, b, c, d

    def get_data_extent(self, reference='center', nodata=None, crs=None):
        """ Return the region characteristics as a tuple (xmin, xmax, ymin,
        ymax). *reference* is a string and may be 'center' or 'edge'. """
        if nodata is None:
            nodata = self.nodata

        if np.isnan(nodata):
            isdata = lambda a: ~np.isnan(a)
        else:
            def isdata(a):
                return a != nodata

        dx, dy = self.transform[2:4]
        sx, sy = self.transform[4:6]
        x0 = x1 = self.transform[0] + 0.5*dx + 0.5*sx
        y0 = y1 = self.transform[1] + 0.5*dy + 0.5*sy
        ny, nx = self.size

        # Reading a row is fast, so process from bottom to top
        bx = x0
        by = y0 + ny*dy
        tx = x0
        ty = y0
        lx = x0 + nx*dx
        ly = y0
        rx = y0
        ry = y0

        jvec = np.arange(nx)

        for i, row in enumerate(self.values):
            x = (x0 + jvec*dx + i*sx)[isdata(row)]
            y = (y0 + i*dy + jvec*sy)[isdata(row)]

            if len(x) != 0:

                xmax, xmin = x.max(), x.min()
                ymax, ymin = y.max(), y.min()
                if ymin < by:
                    by = ymin
                    bx = x[y==ymin]
                if ymax > ty:
                    ty = ymax
                    tx = x[y==ymax]
                if xmin < lx:
                    lx = xmin
                    ly = y[x==xmin]
                if xmax > rx:
                    rx = xmax
                    ry = y[x==xmax]

        if reference == 'center':
            pass
        elif reference == 'edge':
            lx -= 0.5*dx - 0.5*sx
            rx += 0.5*dx + 0.5*sx
            by -= 0.5*dy - 0.5*sy
            ty += 0.5*dy + 0.5*sy
        else:
            raise errors.GridError("`reference` must be 'center' or 'edge'")

        if (crs is not None) and (crs != self.crs):
            lx, ly = crs.project(*self.crs.project(lx, ly, inverse=True))
            rx, ry = crs.project(*self.crs.project(rx, ry, inverse=True))
            bx, by = crs.project(*self.crs.project(bx, by, inverse=True))
            tx, ty = crs.project(*self.crs.project(tx, ty, inverse=True))
        return (lx, rx, by, ty)

    @property
    def data_mask(self):
        """ 8-bit mask of valid data cells """
        if np.isnan(self.nodata):
            isdata = lambda a: ~np.isnan(a)
        else:
            def isdata(a):
                return a != self.nodata
        return isdata(self.values)

    def aschunks(self, size=(-1, -1), overlap=(0, 0), copy=True):
        """ Generator for grid chunks of *size* and *overlap*.
        """
        ny, nx = self.size
        if size == (-1, -1):
            size = (nx//4, ny//4)

        i0 = 0
        j0 = 0

        T0 = self.transform
        while 1:
            if j0 >= nx:
                j0 = 0
                i0 += size[1]-overlap[1]

            if i0 >= ny:
                break

            T = [self.transform[0] + j0*T0[2] + i0*T0[4],
                 self.transform[1] + i0*T0[3] + j0*T0[5],
                 T0[2], T0[3], T0[4], T0[5]]

            yield RegularGrid(T,
                              values=self.values[i0:i0+size[1], j0:j0+size[0]].copy(),
                              crs=self.crs, nodata_value=self.nodata)
            j0 += size[0]-overlap[0]

    def clip(self, xmin, xmax, ymin, ymax, crs=None):
        """ Return a clipped version of grid with cell centers constrained to a
        bounding box defined by *xmin*, *xmax*, *ymin*, *ymax*.

        Optional *crs* argument defines the coordinate reference system of the
        bounding box.
        """
        if crs is not None:
            xg, yg = crs.project([xmin, xmin, xmax, xmax],
                                 [ymin, ymax, ymin, ymax],
                                 inverse=True)
            x, y = self.crs.project(xg, yg)
            xmin = min(x)
            xmax = max(x)
            ymin = min(y)
            ymax = max(y)

        # Convert to center coords
        t = self.transform

        ll = self.get_positions(xmin, ymin)
        lr = self.get_positions(xmax, ymin)
        ul = self.get_positions(xmin, ymax)
        ur = self.get_positions(xmax, ymax)

        i0 = int(np.ceil(min(ll[0], lr[0], ul[0], ur[0])))
        i1 = int(np.floor(max(ll[0], lr[0], ul[0], ur[0]))) + 1
        j0 = int(np.ceil(min(ll[1], lr[1], ul[1], ur[1])))
        j1 = int(np.floor(max(ll[1], lr[1], ul[1], ur[1]))) + 1

        values = self.values[i0:i1,j0:j1].copy()
        x0 = t[0] + j0*t[2] + i0*t[4]
        y0 = t[1] + i0*t[3] + j0*t[5]
        tnew = (x0, y0, t[2], t[3], t[4], t[5])
        return RegularGrid(tnew, values, crs=self.crs, nodata_value=self.nodata)

    def resize(self, bboxnew):
        """ Return a new grid grid with outer edges given by *bboxnew* (xmin,
        ymin, xmax, ymax).
        
        If the grid origin shifts by a non-integer factor of the grid
        resolution, nearest neighbour values are selected. """
        bb = self.bbox
        bbnew = list(bboxnew)
        dx, dy, sx, sy = self.transform[2:]

        # Redefine output bbox so that (xmax-xmin)/dx and (ymax-ymin)/dx are
        # integers
        bbnew[2] = bbnew[0] + dx*math.ceil((bbnew[2]-bbnew[0])/dx)
        bbnew[3] = bbnew[1] + dy*math.ceil((bbnew[3]-bbnew[1])/dy)

        ny, nx = self.size
        nxnew = int((bbnew[2]-bbnew[0])/dx)
        nynew = int((bbnew[3]-bbnew[1])/dy)
        Tnew = [bbnew[0], bbnew[1], dx, dy, sx, sy]
        valnew = self.nodata * np.ones([nynew, nxnew], dtype=self.values.dtype)

        # determine the indices of existing data on the new grid
        j0new = max(0,     int(round((bb[0]-bbnew[0])/dx)))
        j1new = min(nxnew, int(round((bb[2]-bbnew[0])/dx)))
        i0new = max(0,     int(round((bb[1]-bbnew[1])/dy)))
        i1new = min(nynew, int(round((bb[3]-bbnew[1])/dy)))

        # determine the indices of the data to copy on the old grid
        j0 = max(0,  int(round((bbnew[0]-bb[0])/dx)))
        j1 = min(nx, j0+nxnew)
        i0 = max(0,  int(round((bbnew[1]-bb[1])/dy)))
        i1 = min(ny, i0+nynew)

        valnew[i0new:i1new, j0new:j1new] = self.values[i0:i1,j0:j1]
        gridnew = RegularGrid(Tnew, values=valnew, crs=self.crs,
                              nodata_value=self.nodata)
        return gridnew

    def mask_by_poly(self, poly, inplace=False):
        """ Return a grid with all elements outside the bounds of a polygon
        masked by nodata """
        if not poly.isclockwise():
            vertices = poly.get_vertices(self.crs)[::-1]
        else:
            vertices = poly.get_vertices(self.crs)
        x = [a[0] for a in vertices]
        y = [a[1] for a in vertices]

        ny, nx = self.size

        msk = mask_poly(x, y, nx, ny, self.transform)

        if inplace:
            self.values[~msk] = self.nodata
            return self
        else:
            return RegularGrid(self.transform,
                               values=np.where(msk, self.values, self.nodata),
                               crs=self.crs, nodata_value=self.nodata)

    def resample_griddata(self, dx, dy, method='nearest'):
        """ Resample array to have spacing `dx`, `dy' using *scipy.griddata*

        Warning: This method is likely to be removed in the future, once
        the built-in `resample` method matures.

        Parameters
        ----------

        dx : cell dimension, float

        dy : cell dimension, float

        method : interpolation method, string ('nearest', 'linear')
        """
        if not HASSCIPY:
            raise errors.MissingDependencyError("resample_griddata requires scipy")
        ny0, nx0 = self.values.shape[:2]
        dx0, dy0 = self._transform[2:4]
        xllcenter, yllcenter = self.center_llref()
        xurcenter = xllcenter + dx0*nx0
        yurcenter = yllcenter + dy0*ny0
        nx = int((xurcenter - xllcenter) // dx)
        ny = int((yurcenter - yllcenter) // dy)

        xx0, yy0 = self.coordmesh()
        xx, yy = np.meshgrid(np.linspace(xllcenter, xurcenter, nx),
                             np.linspace(yllcenter, yurcenter, ny))
        idata = interpolate.griddata((xx0.flatten(), yy0.flatten()), self.values.flatten(),
                                     (xx.flatten(), yy.flatten()), method=method)
        values = idata.reshape(ny, nx)
        t = self._transform
        tnew = (t[0], t[1], dx, dy, t[4], t[5])
        return RegularGrid(tnew, values=values, crs=self.crs,
                           nodata_value=self.nodata)

    def resample(self, dx, dy, method='nearest'):
        """ Resample array to have spacing `dx`, `dy'. The grid origin remains
        in the same position.

        Parameters
        ----------

        dx : cell dimension, float

        dy : cell dimension, float

        method : interpolation method, string ('nearest',)
        """
        ny0, nx0 = self.values.shape[:2]
        dx0, dy0 = self._transform[2:4]
        xllcenter, yllcenter = self.center_llref()

        if method == 'nearest':
            rx, ry = dx / dx0, dy / dy0
            ny, nx = self.values.shape
            I = np.around(np.arange(ry/2, ny, ry)-0.5).astype(int)
            J = np.around(np.arange(rx/2, nx, rx)-0.5).astype(int)
            if I[-1] == self.values.shape[0]:
                I = I[:-1]
            if J[-1] == self.values.shape[1]:
                J = J[:-1]
            JJ, II = np.meshgrid(J, I)
            values = self.values[II, JJ]
        else:
            raise NotImplementedError('method "{0}" not '
                                      'implemented'.format(method))

        t = self._transform
        tnew = (t[0], t[1], dx, dy, t[4], t[5])
        return RegularGrid(tnew, values=values, crs=self.crs,
                           nodata_value=self.nodata)

    def get_positions(self, x, y):
        """ Return the float column and row indices for the point nearest
        geographical coordinates (x, y). """
        # Calculate this by forming block matrices
        #       | dx sy          |
        #       | sx dy          |
        #   T = |      ....      |
        #       |          dx sy |
        #       |          sx dy |
        #
        #       | x0 + (dx+dy)/2 |
        #       | y0 + (dx+dy)/2 |
        #   S = |    .......     |
        #       | x0 + (dx+dy)/2 |
        #       | y0 + (dx+dy)/2 |
        #
        # where the grid transform is t = (x0, y0, dx, dy, sy, sx)
        #
        # Then the nearest indices J come from solving the system
        # T J = X - S
        try:
            npts = len(x)
            if npts == 1:
                x = x[0]
                y = y[0]
        except TypeError:
            npts = 1

        n = min(70, npts)    # number to process at once
        t = self._transform
        T = (np.diag(np.tile([t[2], t[3]], n)) +
             np.diag(np.tile([t[4], 0], n)[:-1], 1) +
             np.diag(np.tile([t[5], 0], n)[:-1], -1))
        S = np.tile([t[0]+0.5*(t[2]+t[4]), t[1]+0.5*(t[3]+t[5])], n)

        if npts == 1:
            ind = np.linalg.solve(T, np.array([x, y]) - S)
        else:
            ind = np.empty(2*npts, dtype=np.float64)
            i = 0
            while i != npts:
                ip = min(i + n, npts)
                xy = np.zeros(2*n, dtype=np.float64)
                xy[:2*(ip-i)] = np.vstack([x[i:ip], y[i:ip]]).T.ravel()
                ind_ = np.linalg.solve(T, xy - S)
                ind[2*i:2*ip] = ind_[:2*(ip-i)]
                i = ip

        ny, nx = self.size
        j = ind[::2]
        i = ind[1::2]
        return i, j

    def get_indices(self, x, y):
        """ Return the integer column and row indices for the point nearest
        geographical coordinates (`x`, `y`). Raises GridError if points outside
        grid bounding box. """
        ny, nx = self.size
        i, j = self.get_positions(x, y)

        if len(i) != 1:
            i, j = np.round(i).astype(int), np.round(j).astype(int)
            if min(i) < 0 or max(i) > ny-1 or min(j) < 0 or max(j) > nx-1:
                raise errors.GridError("Coordinates outside grid region ({0})".format(self.bbox))
        else:
            i, j = int(round(i[0])), int(round(j[0]))
            if not (0 <= i <= ny-1) or not(0 <= j <= nx-1):
                raise errors.GridError("Coordinates outside grid region ({0})".format(self.bbox))

        return i,j

    def sample_nearest(self, x, y):
        """ Return the value nearest to (`x`, `y`). Nearest grid center
        sampling scheme. """
        i, j = self.get_indices(x, y)
        ny, nx = self.values.shape[:2]
        return self.values[i, j]

    def sample_bilinear(self, x, y):
        """ Return the value nearest to (`x`, `y`). Bilinear sampling scheme.
        """
        i, j = self.get_positions(x, y)
        i0 = np.floor(i).astype(int)
        i1 = np.ceil(i).astype(int)
        j0 = np.floor(j).astype(int)
        j1 = np.ceil(j).astype(int)

        # Handle case where interpolation point is exactly on a grid point
        mska = i0==i1
        mskb = i0==0
        i0[mska&~mskb] -= 1
        i1[mska&mskb] += 1

        mska = j0==j1
        mskb = j0==0
        j0[mska&~mskb] -= 1
        j1[mska&mskb] += 1

        ny, nx = self.size
        if np.any(i1>=ny) or np.any(j1>=nx) or np.any(i0<0) or np.any(j0<0):
            raise errors.GridError("Coordinates outside grid extent({0})"
                    .format(self.get_extent()))

        dx, dy = self._transform[2:4]
        z = (self.values[i0,j0]*(i1-i)*(j1-j) + self.values[i1,j0]*(i-i0)*(j1-j) + \
             self.values[i0,j1]*(i1-i)*(j-j0) + self.values[i1,j1]*(i-i0)*(j-j0))
        return z

    def sample(self, *args, **kwargs):
        """ Return the values nearest positions. Positions may be:

        - a karta.Point instance
        - a karta.Multipoint instance, or
        - a pair of x, y coordinate lists

        Keyword *crs* used when coordinate lists are provided, otherwise the
        coordinate system is assumed to be the same as the grid's coordinate
        system and no transformation is performed.

        Keyword *method* may be one of 'nearest', 'bilinear' (default).
        """
        crs = kwargs.get("crs", None)
        method = kwargs.get("method", "bilinear")

        argerror = TypeError("`grid` takes a Point, a Multipoint, or x, y coordinate lists")
        if hasattr(args[0], "_geotype"):
            crs = args[0].crs
            if args[0]._geotype == "Point":
                x = args[0].x
                y = args[0].y
            elif args[0]._geotype == "Multipoint":
                x, y = args[0].coordinates
            else:
                raise argerror
        else:
            try:
                x = args[0]
                y = args[1]
                if crs is None:
                    crs = self.crs
            except IndexError:
                raise argerror
            if len(x) != len(y):
                raise argerror

        if crs is not self.crs:
            xg, yg = crs.project(x, y, inverse=True)
            x, y = self.crs.project(xg, yg)

        if method == "nearest":
            return self.sample_nearest(x, y)
        elif method == "bilinear":
            return self.sample_bilinear(x, y)
        else:
            raise ValueError("method '{0}' not available".format(method))

    def profile(self, line, resolution=None, **kw):
        """ Sample along a *line* at *resolution*.

        Parameters:
        -----------
        line : `geometry.Line`-like object describing the sampling path

        resolution : sample spacing, taken to be the minimum grid resolution by
        default

        Additional keyword arguments passed to `RegularGrid.sample` (e.g. to
        specify sampling method)

        Returns:
        --------
        vertices : list of (x, y) tuples
        profile : ndarray
        """
        if resolution is None:
            resolution = min(self.transform[2:4])

        remainder = 0
        pt0 = line[0]
        vertices = [pt0.get_vertex(self.crs)]

        for seg in line.segments:
            pos = 0
            az = seg[0].azimuth(seg[1])

            while pos < seg.length:
                distance_to_endpt = pt0.distance(seg[1])
                if distance_to_endpt >= resolution:
                    pt1 = pt0.walk(resolution - remainder, az)
                    pos += resolution - remainder
                    vertices.append(pt1.get_vertex(self.crs))
                    remainder = 0
                    pt0 = pt1
                else:
                    remainder = distance_to_endpt
                    pos = seg.length
                    pt0 = seg[1]

        z = self.sample(*zip(*vertices), **kw)
        return vertices, z

    def as_warpedgrid(self):
        """ Return a copy of grid as a `WarpedGrid`. This is a more general
        grid class that has a larger memory footprint but can represent more
        flexible data layouts.
        """
        Xc, Yc = self.center_coords()
        return WarpedGrid(Xc, Yc, self.values.copy(), crs=self.crs)

    def to_gtiff(self, fnm, compress="PACKBITS", tiled=False, **kw):
        """ Write data to a GeoTiff file using GDAL.
        
        Parameters
        ----------
        fnm: output file name

        compress: "PACKBITS" (default), "DEFLATE", "LZW", "LZMA", or None
        """
        return _gtiff.write(fnm, self, compress=compress, **kw)

    def to_aai(f, reference='corner', nodata_value=-9999):
        """ Save internal data as an ASCII grid. Based on the ESRI standard,
        only isometric grids (i.e. `hdr['dx'] == hdr['dy']` can be saved,
        otherwise `GridIOError` is thrown.

        Parameters:
        -----------
        f : a file-like object or a filename

        reference : specify a header reference ('center' | 'corner')

        nodata_value : specify how NaNs should be represented (int or float)
        """
        if reference not in ('center', 'corner'):
            raise errors.GridIOError("reference in AAIGrid.tofile() must be 'center' or "
                           "'corner'")

        if np.any(self._transform[4:] != 0.0):
            raise errors.GridIOError("ESRI ASCII grids do not support skewed grids")

        ny, nx = self.values.shape[:2]
        x0, y0, dx, dy = self._transform[:4]
        if dx != dy:
            raise errors.GridIOError("ASCII grids require isometric grid cells")

        if not hasattr(f, 'read'):
            f = open(f, "w")

        try:
            data_a = self.values.copy()
            data_a[np.isnan(data_a)] = nodata_value

            f.write("NCOLS {0}\n".format(nx))
            f.write("NROWS {0}\n".format(ny))
            if reference == 'center':
                f.write("XLLCENTER {0}\n".format(x0))
                f.write("YLLCENTER {0}\n".format(y0))
            elif reference == 'corner':
                xllcorner, yllcorner = self.corner_llref()
                f.write("XLLCORNER {0}\n".format(xllcorner))
                f.write("YLLCORNER {0}\n".format(yllcorner))
            f.write("CELLSIZE {0}\n".format(dx))
            f.write("NODATA_VALUE {0}\n".format(nodata_value))
            f.writelines([str(row).replace(',','')[1:-1] +
                            '\n' for row in data_a.tolist()])
        finally:
            f.close()
        return

    def gtiffwrite(self, fnm):
        warnings.warn("method `gtiffwrite` has been renamed `to_gtiff`",
                FutureWarning)
        return self.to_gtiff(fnm)

    def aaiwrite(self, *args, **kwargs):
        warnings.warn("method `aaiwrite` has been renamed `to_aai`",
                FutureWarning)
        return self.to_aai(*args, **kwargs)

class WarpedGrid(Grid):
    """ Warped Grid class. A WarpedGrid contains a fixed number of rows
    and columns and a scalar or vector field defined on the z-axis. Grid
    spacing is not necessarily constant.
    """

    def __init__(self, X, Y, values, crs=None, nodata_value=None):
        """
        Parameters:
        -----------
        X : first-dimension coordinates of grid centers

        Y : second-dimension coordinates of grid centers

        values : dependent m-dimensional quantity (nrows x ncols)
        """

        if any(a is None for a in (X, Y, values)):
            raise errors.GridError('All of (X, Y, values) must be provided')

        if not (X.shape == Y.shape == values.shape[:2]):
            raise errors.GridError('All of (X, Y, values) must share the same '
                            'size over the first two dimensions')

        if crs is None:
            self.crs = Cartesian
        else:
            self.crs = crs

        self.X = X
        self.Y = Y
        self.values = values

        if nodata_value is None:
            self._nodata = get_nodata(self.values.dtype.type)
        else:
            self._nodata = nodata_value
        return

    def __add__(self, other):
        if self._equivalent_structure(other):
            return WarpedGrid(self.X.copy(), self.Y.copy(), self.values+other.values)
        else:
            raise errors.NonEquivalentGridError(self, other)

    def __sub__(self, other):
        if self._equivalent_structure(other):
            return WarpedGrid(self.X.copy(), self.Y.copy(), self.values-other.values)
        else:
            raise errors.NonEquivalentGridError(self, other)

    def _equivalent_structure(self, other):
        return np.all(self.X == other.X) and np.all(self.Y == other.Y) and \
                np.all(self.values.shape == other.values.shape)

    def rotate(self, deg, origin=(0.0, 0.0)):
        """ Rotate grid by *deg* degrees counter-clockwise around *origin*. """
        raise NotImplementedError

    def resample(self, X, Y):
        """ Resample internal grid to the points defined by *X*, *Y*. """
        raise NotImplementedError

def merge(grids, weights=None):
    """ Perform a basic grid merge. Currently limited to grids whose sampling
    is an integer translation from each other.
   
    Parameters:
    -----------

    grids : iterable of Grid objects to combine

    weights : (optional) iterable of weighting factors for computing grid averages
    """

    # Check grid class
    if not all(isinstance(grid, RegularGrid) for grid in grids):
        raise NotImplementedError("All grids must by type RegularGrid")

    T = grids[0].transform
    # Check grid stretch and skew
    for i, grid in enumerate(grids[1:]):
        if grid.transform[2:6] != T[2:6]:
            raise NotImplementedError("grid %d transform stretch and skew "
                    "does not match grid 1" % (i+2,))

    # Check grid offset
    excmsg = "grid %d not an integer translation from grid 1"
    for i, grid in enumerate(grids[1:]):
        if ((grid.transform[0]-T[0]) / float(T[2])) % 1 > 1e-15:
            raise NotImplementedError(excmsg % (i+2,))
        if ((grid.transform[1]-T[1]) / float(T[3])) % 1 > 1e-15:
            raise NotImplementedError(excmsg % (i+2,))

    # Compute weighting coefficients by normalizing *weights*
    if weights is None:
        weights = np.ones(len(grids))
    else:
        weights = np.asarray(weights, dtype=np.float32)

    normalizedweights = weights * len(weights) / weights.sum()

    # Compute final grid extent
    xmin, xmax, ymin, ymax = grids[0].get_extent(reference='edge')
    for grid in grids[1:]:
        _xmin, _xmax, _ymin, _ymax = grid.get_extent(reference='edge')
        xmin = min(xmin, _xmin)
        xmax = max(xmax, _xmax)
        ymin = min(ymin, _ymin)
        ymax = max(ymax, _ymax)

    nx = int(round((xmax-xmin) / T[2]))
    ny = int(round((ymax-ymin) / T[3]))

    # Allocate data array and copy each grid's data
    typ = grids[0].values.dtype.type
    values = np.zeros([ny, nx], dtype=typ)
    counts = np.zeros([ny, nx], dtype=np.float32)
    for grid, w in zip(grids, normalizedweights):
        _xmin, _xmax, _ymin, _ymax = grid.get_extent(reference='edge')
        offx = int((_xmin-xmin) / T[2])
        offy = int((_ymin-ymin) / T[3])
        _ny, _nx = grid.size

        mask = grid.data_mask
        counts[offy:offy+_ny,offx:offx+_nx][mask] += w
        values[offy:offy+_ny,offx:offx+_nx][mask] += typ(grid.values[mask]*w)
        del mask

    validcountmask = (counts!=0.0)
    values[validcountmask] = values[validcountmask] / counts[validcountmask]
    values[~validcountmask] = grids[0].nodata
    Tmerge = [xmin, ymin] + list(T[2:])
    return RegularGrid(Tmerge, values=values, crs=grids[0].crs,
                       nodata_value=grids[0].nodata)

def get_nodata(T):
    """ Return a default value for NODATA given a type

    For unsigned integer types, returns largest representable value.
    For signed integer types, returns smallest negative representable value.
    For floating point types (incl. complex), returns NaN.
    Otherwise, raises ValeError.
    """
    if T in (np.uint8, np.uint16, np.uint32, np.uint64):
        return np.iinfo(T).max
    elif issubclass(T, IntegerType):
        return np.iinfo(T).min
    elif issubclass(T, (numbers.Real, numbers.Complex)):
        return np.nan
    else:
        raise ValueError("No default NODATA value for type {0}".format(T))

def mask_poly(xpoly, ypoly, nx, ny, transform):
    """ Create a grid mask based on a clockwise-oriented polygon.

    Arguments
    ---------
    xpoly: list[float]
    ypoly: list[float]
        sequences of points representing polygon
    nx: int
    ny: int
        size of grid
    transform: list[float]
        affine transformation describing grid layout and origin
        T == [x0, y0, dx, dy, sx, sy]
    """
    mask = np.zeros((ny, nx), dtype=np.int8)

    # find southernmost index (will start and end on this point)
    v = ypoly[0]
    i_bot = 0
    for i in range(1, len(ypoly)):
        if ypoly[i] < v:
            v = ypoly[i]
            i_bot = i

    x0 = xpoly[i_bot]
    y0 = ypoly[i_bot]

    ta, tb, tc, td, te, tf = transform
    i0 = int(round((y0 - min(y0, tb) - tf/tc*(x0 - min(x0, ta))) / (td - tf*te/tc)))
    j0 = int(round((x0 - min(x0, ta) - te/td*(y0 - min(y0, tb))) / (tc - te*tf/td)))

    for el in range(1, len(xpoly)+1):
        idx = (el + i_bot) % len(xpoly)
        x1 = xpoly[idx]
        y1 = ypoly[idx]

        # (Unbounded) grid indices of the segment end points
        i1 = int(round((y1-tb - tf/tc*(x1-ta)) / (td - tf*te/tc)))
        j1 = int(round((x1-ta - te/td*(y1-tb)) / (tc - te*tf/td)))

        # If segment is horizontal or off-grid, ignore
        if ((0 <= i0 < ny) and (0 <= i1 < ny)) or (y1 != y0):

            if y1 > y0:     # mark grid cells to the right

                for i in range(i0, i1):
                    if (0 <= i < ny):
                        j = int(round((i-i0) * (x1-x0)/(y1-y0) + j0))
                        if j < nx:
                            mask[i, max(0, j):] += 1

            else:           # unmark grid cells to the right

                for i in range(i1, i0):
                    if (0 <= i < ny):
                        j = int(round((i-i1) * (x1-x0)/(y1-y0) + j1))
                        if j < nx:
                            mask[i, max(0, j):] -= 1

        x0 = x1
        y0 = y1
        i0 = i1
        j0 = j1

    return mask.astype(np.bool)

