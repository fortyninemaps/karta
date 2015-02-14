""" Classes for basic grid types """

import sys
import copy
import math
import numbers
import numpy as np
from ..crs import Cartesian

IntegerType = (numbers.Integral, np.int32, np.int64)

class Grid(object):
    """ Grid baseclass. Don't use this directly except to implement subclasses.
    The primary attributes defined by all Grid-derived classed are _hdr and
    data.
    """

    @property
    def crs(self):
        return getattr(self, "_crs", None)

    @property
    def nodata(self):
        return self._nodata

    def clipz(self, bounds):
        """ Clip the z-range in place to bounds = [min, max]. """
        self.values = self.values.clip(bounds[0], bounds[1])
        return

    def max(self):
        """ Return the maximum non-nan in self.data. """
        return np.nanmax(self.values)

    def min(self):
        """ Return the minimum non-nan in self.data. """
        return np.nanmin(self.values)

    def minmax(self):
        """ Return the minimum and maximum value of data array. """
        return (self.min(), self.max())

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

    in which *(a,b)* defines the grid origin, *(c,d)* defines the resolution in
    the horizontal and vertical directions, and *(e,f)* can be used to define a
    rotation. In the common case of an orthogonal grid with north "up",

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
            self._transform = tuple([transform[f] for f in
                    ("xllcenter", "yllcenter", "dx", "dy", "xrot", "yrot")])
        elif len(transform) == 6:
            self._transform = tuple(transform)
        elif len(transform) == 4:
            self._transform = (transform[0], transform[1],
                               transform[2], transform[3], 0, 0)
        else:
            raise GridError("RegularGrid must be initialized with a transform"
                            "iterable or dictionary")

        if values is not None:
            self.values = values
        else:
            self.values = np.atleast_2d([np.nan])

        if crs is None:
            self._crs = Cartesian
        else:
            self._crs = crs

        if nodata_value is None:
            self._nodata = get_nodata(self.values.dtype.type)
        else:
            self._nodata = nodata_value
        return

    def __add__(self, other):
        if self._equivalent_structure(other):
            return RegularGrid(copy.copy(self.transform),
                               values=self.values+other.values)
        else:
            raise NonEquivalentGridError(self, other)

    def __sub__(self, other):
        if self._equivalent_structure(other):
            return RegularGrid(copy.copy(self.transform),
                               values=self.values-other.values)
        else:
            raise NonEquivalentGridError(self, other)

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
        return self._transform[0], self._transform[1]

    def corner_llref(self):
        """ Return the 'lower-left' reference in terms of a center coordinate.
        """
        xllcorner = self._transform[0] - 0.5 * self._transform[2]
        yllcorner = self._transform[1] - 0.5 * self._transform[3]
        return (xllcorner, yllcorner)

    def coordmesh(self, *args, **kwargs):
        """ Alias for center_coords() """
        return self.center_coords(*args, **kwargs)

    def center_coords(self):
        t = self._transform
        xcoords = np.empty(self.values.shape[:2])
        ycoords = np.empty(self.values.shape[:2])
        irow = np.arange(self.values.shape[0])
        for i in range(self.values.shape[1]):
            xcoords[:,i] = t[0] + i*t[2] + irow*t[4]
            ycoords[:,i] = (t[1] + irow*t[3] + i*t[5])
        return xcoords, ycoords

    def vertex_coords(self):
        """ Return the coordinates of vertices. """
        xllcorner, yllcorner = self.corner_llref()
        t = self._transform
        ny, nx = self.values.shape[:2]
        xurcorner = xllcorner + (nx+1) * t[2] + (ny+1) * t[4]
        yurcorner = yllcorner + (ny+1) * t[3] + (nx+1) * t[5]
        return np.meshgrid(np.linspace(xllcorner, xurcorner, nx + 1),
                           np.linspace(yurcorner, yllcorner, ny + 1))

    @property
    def bbox(self):
        extents = self.get_extents(reference="edge")
        return extents[0], extents[2], extents[1], extents[3]

    def get_extents(self, reference='center'):
        """ Return the region characteristics as a tuple (xmin, xmax, ymin,
        ymax). *reference* may be `center` or `edge`. """
        if reference == 'center':
            x0, y0 = self.center_llref()
            n = -1
        elif reference == 'edge':
            x0, y0 = self.corner_llref()
            n = 0
        else:
            raise GridError("`reference` must be 'center' or 'edge'")
        ny, nx = self.values.shape[:2]
        dx, dy = self._transform[2:4]
        sx, sy = self._transform[4:]

        sgn = lambda a: 0 if a==0 else a/abs(a)
        if sgn(dx) == sgn(sx):
            xrng = dx*(nx+n) + sx*(ny+n)
        else:
            xrng = dx*(nx+n) - sx*(ny+n)

        if sgn(dy) == sgn(sy):
            yrng = dy*(ny+n) + sy*(nx+n)
        else:
            yrng = dy*(ny+n) - sy*(nx+n)
        return min(x0, x0+xrng), max(x0, x0+xrng), min(y0, y0+yrng), max(y0, y0+yrng)

    def clip(self, xmin, xmax, ymin, ymax, crs=None):
        """ Return a clipped version of grid constrained to a bounding box. """
        if crs is not None:
            xg, yg = crs.proj([xmin, xmax], [ymin, ymax], inverse=True)
            x, y = self.crs.proj(xg, yg)
            xmin, xmax = x
            ymin, ymax = y
        else:
            crs = self.crs

        ll = self.get_indices(xmin, ymin)
        lr = self.get_indices(xmax, ymin)
        ul = self.get_indices(xmin, ymax)
        ur = self.get_indices(xmax, ymax)

        i0 = min(ll[0], lr[0], ul[0], ur[0])
        i1 = max(ll[0], lr[0], ul[0], ur[0])
        j0 = min(ll[1], lr[1], ul[1], ur[1])
        j1 = max(ll[1], lr[1], ul[1], ur[1])

        values = self.values[i0:i1,j0:j1].copy()
        t = self.transform
        x0 = t[0] + j0*t[2] + i0*t[4]
        y0 = t[1] + i0*t[3] + j0*t[5]
        tnew = (x0, y0, t[2], t[3], t[4], t[5])
        return RegularGrid(tnew, values, crs=self.crs)

    def resample_griddata(self, dx, dy, method='nearest'):
        """ Resample array to have spacing `dx`, `dy' using *scipy.griddata*

        Parameters
        ----------

        dx : cell dimension, float

        dy : cell dimension, float

        method : interpolation method, string ('nearest', 'linear')
        """
        from scipy.interpolate import griddata
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
        idata = griddata((xx0.flatten(), yy0.flatten()), self.values.flatten(),
                         (xx.flatten(), yy.flatten()), method=method)
        values = idata.reshape(ny, nx)
        t = self._transform
        tnew = (t[0], t[1], dx, dy, t[4], t[5])
        return RegularGrid(tnew, values)

    def resample(self, dx, dy, method='nearest'):
        """ Resample array to have spacing `dx`, `dy'.

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
            I = np.around(np.arange(ry/2, self.values.shape[0], ry)).astype(int)
            J = np.around(np.arange(rx/2, self.values.shape[1], rx)).astype(int)
            JJ, II = np.meshgrid(J, I)
            values = self.values[II, JJ]
        else:
            raise NotImplementedError('method "{0}" not '
                                      'implemented'.format(method))

        t = self._transform
        tnew = (t[0], t[1], dx, dy, t[4], t[5])
        return RegularGrid(tnew, values)

    def get_positions(self, x, y):
        """ Return the column and row indices for the point nearest
        geographical coordinates (x, y). """
        # Calculate this by forming block matrices
        #       | dx sy          |
        #       | sx dy          |
        #   T = |      ....      |
        #       |          dx sy |
        #       |          sx dy |
        #
        #       | x0 |
        #       | y0 |
        #   S = | .. |
        #       | x0 |
        #       | y0 |
        #
        # where the grid transform is t = (x0, y0, dx, dy, sy, sx)
        #
        # Then the nearest indices J come from solving the system
        # T J = X - S
        try:
            npts = len(x)
        except TypeError:
            npts = 1

        n = min(70, npts)    # number to process at once
        t = self._transform
        T = (np.diag(np.tile([t[2], t[3]], n)) +
             np.diag(np.tile([t[4], 0], n)[:-1], 1) +
             np.diag(np.tile([t[5], 0], n)[:-1], -1))
        S = np.tile([t[0], t[1]], n)

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
        #j = ind[::2].clip(0, nx-1)
        #i = ind[1::2].clip(0, ny-1)
        j = ind[::2]
        i = ind[1::2]
        return i, j

    def get_indices(self, x, y):
        """ Return the column and row indices for the point nearest
        geographical coordinates (x, y). """
        ny, nx = self.size
        i, j = self.get_positions(x, y)

        if len(i) != 1:
            i, j = np.round(i).astype(int), np.round(j).astype(int)
            if min(i) < 0 or max(i) > ny-1 or min(j) < 0 or max(j) > nx-1:
                raise GridError("Coordinates outside grid region ({0})".format(self.bbox))
        else:
            i, j = int(round(i[0])), int(round(j[0]))
            if not (0 <= i <= ny-1) or not(0 <= j <= nx-1):
                raise GridError("Coordinates outside grid region ({0})".format(self.bbox))

        return i,j

    def sample_nearest(self, x, y):
        """ Return the value nearest to (`x`, `y`). Nearest grid center
        sampling scheme. """
        i, j = self.get_indices(x, y)
        ny, nx = self.values.shape[:2]
        return self.values[i, j]

    def sample_bilinear(self, x, y):
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

        dx, dy = self._transform[2:4]
        z = (self.values[i0,j0]*(i1-i)*(j1-j) + self.values[i1,j0]*(i-i0)*(j1-j) + \
             self.values[i0,j1]*(i1-i)*(j-j0) + self.values[i1,j1]*(i-i0)*(j-j0))
        return z

    def sample(self, x, y, crs=None, method="bilinear"):
        """ Return the values nearest (`x`, `y`), where `x` and `y` may be
        equal length vectors. *method* may be one of `nearest`, `bilinear`. """
        if crs is not None:
            xg, yg = crs.project(x, y, inverse=True)
            x, y = self.crs.project(xg, yg)

        if method == "nearest":
            return self.sample_nearest(x, y)
        elif method == "bilinear":
            return self.sample_bilinear(x, y)
        else:
            raise ValueError("method \"{0}\" not available".format(method))

    def profile(self, line, resolution=None, **kw):
        """ Sample along a line defined as `segments`.

        Parameters:
        -----------
        line : Line-like object describing the sampling path
        resolution : sample spacing

        Additional keyword arguments passed to `RegularGrid.sample`

        Returns:
        --------
        profile : ndarray
        """
        if resolution is None:
            resolution = min(self.transform[2:4])

        np = math.ceil(line.length / resolution)
        vertices = []

        remainder = 0
        pt0 = line[0]
        vertices.append(pt0.get_vertex(self.crs))

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
        """ Return a copy as a WarpedGrid instance. This is a more general
        grid class that has a larger memory footprint but can represent more
        flexible data layouts. """
        Xc, Yc = self.center_coords()
        return WarpedGrid(Xc, Yc, self.values.copy())

    def aaiwrite(self, f, reference='corner', nodata_value=-9999):
        """ Save internal data as an ASCII grid. Based on the ESRI standard,
        only isometric grids (i.e. `hdr['dx'] == hdr['dy']` can be saved,
        otherwise `GridIOError` is thrown.

        Parameters:
        -----------
        f : either a file-like object or a filename
        reference : specify a header reference ("center" | "corner")
        nodata_value : specify how nans should be represented (int or float)
        """
        if reference not in ('center', 'corner'):
            raise GridIOError("reference in AAIGrid.tofile() must be 'center' or "
                           "'corner'")

        if np.any(self._transform[4:] != 0.0):
            raise GridIOError("ESRI ASCII grids do not support skewed grids")

        ny, nx = self.values.shape[:2]
        x0, y0, dx, dy = self._transform[:4]
        if dx != dy:
            raise GridIOError("ASCII grids require isometric grid cells")

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

class WarpedGrid(Grid):
    """ Warped Grid class. A WarpedGrid contains a fixed number of rows
    and columns and a scalar or vector field defined on the z-axis. Grid
    spacing is not necessarily constant.
    """

    def __init__(self, X, Y, values, nodata_value=None):
        """
        Parameters:
        -----------
        X : first-dimension coordinates of grid centers
        Y : second-dimension coordinates of grid centers
        values : dependent m-dimensional quantity (nrows x ncols)
        """

        if any(a is None for a in (X, Y, values)):
            raise GridError('All of (X, Y, values) must be provided')

        if not (X.shape == Y.shape == values.shape[:2]):
            raise GridError('All of (X, Y, values) must share the same '
                            'size over the first two dimensions')

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
            raise NonEquivalentGridError(self, other)

    def __sub__(self, other):
        if self._equivalent_structure(other):
            return WarpedGrid(self.X.copy(), self.Y.copy(), self.values-other.values)
        else:
            raise NonEquivalentGridError(self, other)

    def _equivalent_structure(self, other):
        return np.all(self.X == other.X) and np.all(self.Y == other.Y) and \
                np.all(self.values.shape == other.values.shape)

    def rotate(self, deg, origin=(0.0, 0.0)):
        """ Rotate grid by *deg* degrees counter-clockwise around *origin*. """
        raise NotImplementedError

    def resample(self, X, Y):
        """ Resample internal grid to the points defined by `X`, `Y`. """
        raise NotImplementedError


class GridError(Exception):
    def __init__(self, message=''):
        self.message = message
    def __str__(self):
        return self.message


class GridIOError(GridError):
    def __init__(self, message=''):
        self.message = message


class NonEquivalentGridError(GridError):
    def __init__(self, A, B, message=''):
        if len(message) == 0:
            self.message = ("{0} and {1} do not share equivalent "
                            "grid layouts".format(A, B))
        else:
            self.message = message

def get_nodata(T):
    """ Return a default value for NODATA given a type (e.g. int, float,
    complex).
    """
    if issubclass(T, IntegerType):
        return -9999
    elif issubclass(T, (numbers.Real, numbers.Complex)):
        return np.nan
    else:
        raise ValueError("No default NODATA value for type {0}".format(T))

