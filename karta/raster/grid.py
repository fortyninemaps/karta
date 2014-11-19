""" Classes for basic grid types """

import sys
import copy
from math import sqrt
import numpy as np
from . import _aai          # Contains the ascread driver
from ..crs import crsreg, CRS, pyproj

try:
    from . import _gtiff
    HAS_GDAL = True
except ImportError:
    HAS_GDAL = False

try:
    from .fill_sinks import fill_sinks
except ImportError:
    sys.stderr.write("Compiled fill_sinks not available. Falling back to "
                     "Python version\n")
    from .raster import fill_sinks

class Grid(object):
    """ Grid baseclass. Don't use this directly except to implement subclasses.
    The primary attributes defined by all Grid-derived classed are _hdr and
    data.
    """

    @property
    def crs(self):
        return getattr(self, "_crs", None)

    def clipz(self, bounds):
        """ Clip the z-range in place to bounds = [min, max]. """
        self.Z = self.Z.clip(bounds[0], bounds[1])
        return

    def max(self):
        """ Return the maximum non-nan in self.data. """
        return np.nanmax(self.Z)

    def min(self):
        """ Return the minimum non-nan in self.data. """
        return np.nanmin(self.Z)

    def minmax(self):
        """ Return the minimum and maximum value of data array. """
        return (self.min(), self.max())

    def copy(self):
        """ Return a copy. """
        return copy.deepcopy(self)


class RegularGrid(Grid):
    """ Regular (structured) grid class. A RegularGrid contains a fixed number
    of rows and columns with a constant spacing and a scalar or vector field
    defined as `Z`.

    Positions on a RegularGrid are referenced using an affine transform such
    that
    
    X = a + j * c + i * e
    Y = b + i * d + j * f

    in which (a,b) define the grid origin, (c,d) define the resolution in the
    horizontal and vertical directions, and (e,f) can be used to define a
    rotation. In the usual case of an orthogonal grid with north "up", e = f =
    0.
    """
    def __init__(self, transform, Z=None, crs=None):
        """ Create a RegularGrid instance with cells referenced according to
        *transform*, which is an iterable or dictionary consisting of
        (xllcenter, yllcenter, dx, dy, xrot, yrot).

        Optional parameters:
        --------------------
        Z : grid data (nrows x ncols) [default NaN]
        crs : Grid coordinate reference system [default None]
        """
        if hasattr(transform, "keys"):
            self._transform = tuple([transform[f] for f in
                            ("xllcenter", "yllcenter", "dx", "dy", "xskew", "yskew")])
        elif len(transform) == 6:
            self._transform = tuple(transform)
        elif len(transform) == 4:
            self._transform = (transform[0], transform[1],
                               transform[2], transform[3], 0, 0)
        else:
            raise GridError("RegularGrid must be initialized with a transform"
                            "iterable or dictionary")

        if Z is not None:
            self.Z = Z
        else:
            self.Z = np.atleast_2d([np.nan])

        if crs is None:
            self._crs = crsreg.CARTESIAN
        else:
            self._crs = crs
        return

    def __add__(self, other):
        if self._equivalent_structure(other):
            return RegularGrid(copy.copy(self.transform),
                               Z=self.Z+other.Z)
        else:
            raise NonEquivalentGridError(self, other)

    def __sub__(self, other):
        if self._equivalent_structure(other):
            return RegularGrid(copy.copy(self.transform),
                               Z=self.Z-other.Z)
        else:
            raise NonEquivalentGridError(self, other)

    def _equivalent_structure(self, other):
        return (self._transform == other._transform) and \
               (self.Z.shape == other.Z.shape)

    @property
    def transform(self):
        return self._transform

    @property
    def size(self):
        return self.Z.shape

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
        xcoords = np.empty(self.Z.shape[:2])
        ycoords = np.empty(self.Z.shape[:2])
        irow = np.arange(self.Z.shape[0])
        for i in range(self.Z.shape[1]):
            xcoords[:,i] = t[0] + i*t[2] + irow*t[4]
            ycoords[:,i] = (t[1] + irow*t[3] + i*t[5])[::-1]
        return xcoords, ycoords

    # This version is too space-ineffficient, and results in MemoryErrors on
    # moderately-sized grids
    # def center_coords(self):
    #     """ Return the cell-center coordinates. """
    #     ii, jj = np.meshgrid(np.arange(self.Z.shape[1]),
    #                          np.arange(self.Z.shape[0]))
    #     t = self._transform
    #     return (t[0] + ii * t[2] + jj[::-1] * t[4],
    #             t[1] + jj[::-1] * t[3] + ii * t[5])

    def vertex_coords(self):
        """ Return the coordinates of vertices. """
        xllcorner, yllcorner = self.corner_llref()
        t = self._transform
        ny, nx = self.Z.shape[:2]
        xurcorner = xllcorner + (nx+1) * t[2] + (ny+1) * t[4]
        yurcorner = yllcorner + (ny+1) * t[3] + (nx+1) * t[5]
        return np.meshgrid(np.linspace(xllcorner, xurcorner, nx + 1),
                           np.linspace(yurcorner, yllcorner, ny + 1))

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
        ny, nx = self.Z.shape[:2]
        dx, dy = self._transform[2:4]
        return (x0, x0 + dx*(nx+n), y0, y0 + dy*(ny+n))

    def resample_griddata(self, dx, dy, method='nearest'):
        """ Resample array to have spacing `dx`, `dy' using *scipy.griddata*

        Parameters
        ----------

        dx : cell dimension, float

        dy : cell dimension, float

        method : interpolation method, string ('nearest', 'linear')
        """
        from scipy.interpolate import griddata
        ny0, nx0 = self.Z.shape[:2]
        dx0, dy0 = self._transform[2:4]
        xllcenter, yllcenter = self.center_llref()
        xurcenter = xllcenter + dx0*nx0
        yurcenter = yllcenter + dy0*ny0
        nx = int((xurcenter - xllcenter) // dx)
        ny = int((yurcenter - yllcenter) // dy)

        xx0, yy0 = self.coordmesh()
        xx, yy = np.meshgrid(np.linspace(xllcenter, xurcenter, nx),
                             np.linspace(yllcenter, yurcenter, ny))
        idata = griddata((xx0.flatten(), yy0.flatten()), self.Z.flatten(),
                         (xx.flatten(), yy.flatten()), method=method)
        Z = idata.reshape(ny, nx)[::-1]
        t = self._transform
        tnew = (t[0], t[1], dx, dy, t[4], t[5])
        return RegularGrid(tnew, Z)

    def resample(self, dx, dy, method='nearest'):
        """ Resample array to have spacing `dx`, `dy'.

        Parameters
        ----------

        dx : cell dimension, float

        dy : cell dimension, float

        method : interpolation method, string ('nearest',)
        """
        ny0, nx0 = self.Z.shape[:2]
        dx0, dy0 = self._transform[2:4]
        xllcenter, yllcenter = self.center_llref()

        if method == 'nearest':
            rx, ry = dx / dx0, dy / dy0
            I = np.around(np.arange(ry/2, self.Z.shape[0], ry)).astype(int)
            J = np.around(np.arange(rx/2, self.Z.shape[1], rx)).astype(int)
            JJ, II = np.meshgrid(J, I)
            Z = self.Z[II, JJ]
        else:
            raise NotImplementedError('method "{0}" not '
                                      'implemented'.format(method))

        t = self._transform
        tnew = (t[0], t[1], dx, dy, t[4], t[5])
        return RegularGrid(tnew, Z)

    # def resize(self, te):
    #     """ Resize array to fit within extents given by te. If the new
    #     dimensions are smaller, the data is clipped. If they are larger,
    #     nan padding is added.

    #     *te*        :   tuple of center coordinates in the form
    #                     (xmin, xmax, ymin, ymax).

    #     Returns None.
    #     """
    #     t = self._transform
    #     ny, nx = self.Z.shape[:2]
    #     xll1, yll1 = self.center_llref()
    #     xur1 = xll1 + t[2] * (nx - 1)
    #     yur1 = yll1 + t[3] * (ny - 1)

    #     xmin2 = te[0]
    #     xmax2 = te[1]
    #     ymin2 = te[2]
    #     ymax2 = te[3]

    #     data_a = self.data.copy()
    #     ny = data_a.shape[0]

    #     # The left side
    #     Dx = int(np.floor((xmin2-xmin1) / float(self._hdr['dx'])))
    #     if Dx < 0:
    #         data_a = np.hstack([np.nan*np.ones([ny, Dx]), data_a])
    #     elif Dx > 0:
    #         data_a = data_a[:,Dx:]
    #     xllcenter = xmin1 + self._hdr['dx'] * Dx

    #     # The right side
    #     Dx = int(np.floor((xmax2-xmax1) / float(self._hdr['dx'])))
    #     if Dx > 0:
    #         data_a = np.hstack([data_a, np.nan*np.ones([ny, Dx])])
    #     elif Dx < 0:
    #         data_a = data_a[:,:Dx]
    #     self._hdr['nx'] = data_a.shape[1]

    #     _, nx = data_a.shape

    #     # The bottom
    #     Dy = int(np.ceil((ymin2-ymin1) / float(self._hdr['dy'])))
    #     if Dy < 0:
    #         data_a = np.vstack([data_a, np.nan*np.ones([Dy, nx])])
    #     elif Dy > 0:
    #         data_a = data_a[Dy:,:]
    #     yllcenter = ymin1 + self._hdr['dy'] * Dy

    #     # The top
    #     Dy = int(np.floor((ymax2-ymax1) / float(self._hdr['dy'])))
    #     if Dy > 0:
    #         data_a = np.vstack([np.nan*np.ones([Dy, nx]), data_a])
    #     elif Dy < 0:
    #         data_a = data_a[:Dy,:]
    #     self._hdr['ny'] = data_a.shape[0]

    #     self.data = data_a
    #     self._hdr['xllcenter'] = xllcenter
    #     self._hdr['yllcenter'] = yllcenter

    #     return

    def get_indices(self, x, y):
        """ Return the column and row indices for the point nearest
        geographical coordinates (x, y). """
        # Calculate this by forming matrices
        #       | dx sy |
        #   T = |       |
        #       | sx dy |
        #
        #       | x0 |
        #   S = |    |
        #       | y0 |
        #
        # where the grid transform is t = (x0, y0, dx, dy, sy, sx)
        #
        # Then the nearest indices J come from solving the system
        # T J = X - S
        ny, nx = self.size
        t = self._transform
        T = np.array([[t[2], t[4]], [t[5], t[3]]])
        S = np.array([t[0], t[1]])
        ind = np.linalg.solve(T, np.array([x, y])-S)
        ind = ind.clip((0, 0), (nx-1, ny-1))
        return np.round(ind).astype(int)

    def sample_nearest(self, x, y):
        """ Return the value nearest to (`x`, `y`). Nearest grid center
        sampling scheme. """
        xi, yi = self.get_indices(x, y)
        ny, nx = self.Z.shape[:2]
        if (np.any(xi < 0) or np.any(xi >= nx) or
            np.any(yi < 0) or np.any(yi >= ny)):
            raise GridError("coordinates ({0}, {1}) are outside grid region "
                            "({2}, {3})".format(xi, yi, nx, ny))
        return self.Z[yi, xi]

    def sample(self, x, y, method="nearest"):
        """ Return the values nearest (`x`, `y`), where `x` and `y` may be
        equal length vectors. *method* may be one of `nearest`, `linear`. """
        if method == "nearest":
            return self.sample_nearest(x, y)
        elif method == "linear":
            from scipy.interpolate import griddata
            Xd, Yd = self.center_coords()
            return griddata((Xd.flat, Yd.flat), self.Z.flat,
                            (x,y), method="linear")
        else:
            raise ValueError("method \"{0}\" not understood".format(method))

    def get_profile(self, segments, resolution=10.0):
        """ Sample along a line defined as `segments`. Does not interpolate.

        Parameters:
        -----------
        segments : iterable containing (x,y) pairs
        resolution : sample spacing

        Returns:
        --------
        profile : ndarray
        """

        z = []
        p = 0
        for s, f in zip(segments[:-1], segments[1:]):

            x0 = s[0]
            y0 = s[1]
            xf = f[0]
            yf = f[1]

            xlen = xf-x0
            ylen = yf-y0
            d = sqrt((xf - x0)**2 + (yf - y0)**2)

            while p < d:
                fd = p / d
                fx = fd*xlen
                fy = fd*ylen
                xi, yi = self.get_indices(x0+fx, y0+fy)
                z.append(self.Z[yi, xi])
                p += resolution
            p -= d

        return np.array(z)

    def fill_sinks(self):
        """ Fill depressions. Use the algorithm of Wang and Liu (2006). """
        return RegularGrid(self._transform, Z=fill_sinks(self.Z))

    def as_warpedgrid(self):
        """ Return a copy as a WarpedGrid instance. This is a more general
        grid class that has a larger memory footprint but can represent more
        flexible data layouts. """
        Xc, Yc = self.center_coords()
        return WarpedGrid(Xc, Yc, self.Z.copy())

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

        ny, nx = self.Z.shape[:2]
        x0, y0, dx, dy = self._transform[:4]
        if dx != dy:
            raise GridIOError("ASCII grids require isometric grid cells")

        if not hasattr(f, 'read'):
            f = open(f, "w")

        try:
            data_a = self.Z.copy()
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

    def __init__(self, X, Y, Z):
        """
        Parameters:
        -----------
        X : first-dimension coordinates of grid centers
        Y : second-dimension coordinates of grid centers
        Z : dependent m-dimensional quantity (nrows x ncols)
        """

        if any(a is None for a in (X, Y, Z)):
            raise GridError('All of (X, Y, Z) must be provided')

        if not (X.shape == Y.shape == Z.shape[:2]):
            raise GridError('All of (X, Y, Z) must share the same '
                            'size over the first two dimensions')

        self.X = X
        self.Y = Y
        self.Z = Z
        return

    def __add__(self, other):
        if self._equivalent_structure(other):
            return WarpedGrid(self.X.copy(), self.Y.copy(), self.Z+other.Z)
        else:
            raise NonEquivalentGridError(self, other)

    def __sub__(self, other):
        if self._equivalent_structure(other):
            return WarpedGrid(self.X.copy(), self.Y.copy(), self.Z-other.Z)
        else:
            raise NonEquivalentGridError(self, other)

    def _equivalent_structure(self, other):
        return np.all(self.X == other.X) and np.all(self.Y == other.Y) and \
                np.all(self.Z.shape == other.Z.shape)

    def rotate(self, deg, origin=(0.0, 0.0)):
        """ Rotate grid by *deg* degrees counter-clockwise around *origin*. """
        raise NotImplementedError

    def resample(self, X, Y):
        """ Resample internal grid to the points defined by `X`, `Y`. """
        raise NotImplementedError

    def fill_sinks(self):
        """ Fill depressions. Use the algorithm of Wang and Liu (2006). """
        return WarpedGrid(self.X.copy(), self.Y.copy(), fill_sinks(self.Z))


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

def aairead(fnm):
    """ Convenience function to open a ESRI ASCII grid and return a RegularGrid
    instance.
    """
    Z, aschdr = _aai.aairead(fnm)
    t = {'xllcenter': aschdr['xllcorner'] + 0.5 * aschdr['cellsize'],
         'yllcenter': aschdr['yllcorner'] + 0.5 * aschdr['cellsize'],
         'dx'       : aschdr['cellsize'],
         'dy'       : aschdr['cellsize'],
         'xskew'    : 0.0,
         'yskew'    : 0.0}
    Z[Z==aschdr['nodata_value']] = np.nan
    return RegularGrid(t, Z=Z)

def gtiffread(fnm):
    """ Convenience function to open a GeoTIFF and return a RegularGrid
    instance.
    """
    if not HAS_GDAL:
        raise NotImplementedError("Right now, loading GeoTiffs requires GDAL.")
    arr, hdr = _gtiff.read(fnm)
    t = {'xllcenter'  : hdr['xulcorner'] + 0.5 * (hdr['dx'] + hdr['sx']),
         'yllcenter'  : hdr['yulcorner'] + (hdr['ny'] - 0.5) * hdr['dy'] - 0.5 * hdr['sy'],
         'dx'         : hdr['dx'],
         'dy'         : -hdr['dy'],
         'xskew'      : hdr['sx'],
         'yskew'      : hdr['sy']}

    if "lonlat" in hdr["srs"]["proj4"]:
        crstype = "geographical"
    else:
        crstype = "projected"

    crs = CRS(proj=pyproj.Proj(hdr["srs"]["proj4"]),
              geod=pyproj.Geod(a=hdr["srs"]["semimajor"],
                               f=hdr["srs"]["flattening"]),
              crstype=crstype)
    return RegularGrid(t, Z=arr.squeeze(), crs=crs)

