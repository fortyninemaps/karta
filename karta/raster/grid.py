"""
Raster grid representations
"""
import copy
import math
import numbers
import warnings
import numpy as np
from . import _gdal
from . import crfuncs
from .band import SimpleBand, CompressedBand, BandIndexer
from .coordgen import CoordinateGenerator
from .. import errors
from ..crs import Cartesian

BAND_CLASS_DEFAULT = CompressedBand
CRS_DEFAULT = Cartesian

class Grid(object):
    """ Grid base class """

    @property
    def nodata(self):
        return self._nodata

    def max(self):
        """ Return the maximum non-nan in self.data """
        tmp = self[self.data_mask]
        if len(tmp) != 0:
            return tmp.max()
        else:
            return np.nan

    def min(self):
        """ Return the minimum non-nan in self.data """
        tmp = self[self.data_mask]
        if len(tmp) != 0:
            return tmp.min()
        else:
            return np.nan

    def minmax(self):
        """ Return the minimum and maximum value of data array """
        tmp = self[self.data_mask]
        if len(tmp) != 0:
            return (tmp.min(), tmp.max())
        else:
            return (np.nan, np.nan)

    def apply(self, func, copy=False):
        """ Apply a function *func* to grid values """
        msk = self.data_mask_full
        if copy:
            newgrid = self.copy()
            newgrid[:,:] = np.where(msk, func(self[:,:]), self._nodata)
            return newgrid
        else:
            self[:,:] = np.where(msk, func(self[:,:]), self._nodata)
            return self

    def copy(self):
        """ Return a deep copy """
        return copy.deepcopy(self)


class RegularGrid(Grid):
    """ Regular  structured grid class. A RegularGrid contains a fixed number
    of rows and columns with a constant spacing and a list of bands
    representing scalar fields.

    Positions on a RegularGrid are referenced to their array indices (i,j)
    using an affine transform *T* such that

    ::

        T = (a, b, c, d, e, f)

        x = a + j*c + i*e
        y = b + i*d + j*f

    or as a matrix transformation,

    ::

        |e  c  a|   |i|   |x|
        |       | * |j| = | |
        |d  f  b|   |1|   |y|

    where (x,y) are the coordinates of the lower left corner of each pixel.

    Meaning:

    - (a,b) defines the lower-left grid cell corner
    - (c,d) defines the resolution in the horizontal and vertical directions
    - (e,f) can be used to define a rotation

    In the common case of a "east-right, north-up" grid, e = f = 0.
    """
    def __init__(self, transform, values=None, bands=None, crs=None,
            nodata_value=None, bandclass=None):
        """ Create a RegularGrid instance.

        Parameters
        ----------
        transform : 6-tuple of floats
            a six element list, tuple, or dictionary defining the geotransform
            of the grid: ``[xllcorner, yllcorner, dx, dy, sx, sy]``
        values : ndarray (nrows x ncols), optional
            grid data from which raster bands will be formed. A two-dimensional
            array will be transformed into a single band, while a
            three-dimension MxNxP array will become P MxN bands. The structure
            used to store band data internally can be controlled by the
            *bandclass* keyword argument.
        bands : list of band types, optional
            list of allocated Band-like objects. If provided, they are assumed
            to be of a single type, which overrides *bandclass*.
        crs : karta.crs.CRS subclass, optional
            Grid coordinate reference system. default CRS_DEFAULT
        nodata_value : number, optional
            specifies a value to represent NoData or null pixels. The default
            is chosen depending on the input type of *values* or *bands*. If
            neither is provided, the default is NaN.
        bandclass : class, optional
            indicates the band class used to represent grid data. default
            BAND_CLASS_DEFAULT
        """

        if hasattr(transform, "keys"):
            self._transform = tuple([float(transform[f]) for f in
                    ("xllcorner", "yllcorner", "dx", "dy", "sx", "sy")])
        elif len(transform) == 6:
            self._transform = tuple(float(a) for a in transform)
        elif len(transform) == 4:
            self._transform = (float(transform[0]), float(transform[1]),
                               float(transform[2]), float(transform[3]), 0, 0)
        else:
            raise errors.GridError("RegularGrid must be initialized with a "
                                   " transform iterable or dictionary")

        if bands is not None:
            self._bndcls = type(bands[0])
        elif bandclass is None:
            self._bndcls = BAND_CLASS_DEFAULT
        else:
            self._bndcls = bandclass

        if bands is not None:
            self.bands = bands
        else:
            self.bands = []

        if bands is None and (values is not None):
            if values.ndim == 2:
                band = self._bndcls(values.shape, values.dtype.type)
                band[:,:] = values
                self.bands.append(band)
            elif values.ndim == 3:
                for ibnd in range(values.shape[2]):
                    band = self._bndcls(values.shape[:2], values.dtype.type)
                    band[:,:] = values[:,:,ibnd]
                    self.bands.append(band)
            else:
                raise ValueError("`values` must have two or three dimensions")

        self._bandindexer = BandIndexer(self.bands)

        if crs is None:
            self.crs = CRS_DEFAULT
        else:
            self.crs = crs

        if nodata_value is None:
            if len(self.bands) == 0:
                self._nodata = np.nan
            else:
                self._nodata = get_nodata(self.bands[0].dtype)
        else:
            self._nodata = nodata_value
        return

    def __add__(self, other):
        if self._equivalent_structure(other):
            return RegularGrid(copy.copy(self.transform),
                               values=self[:,:]+other[:,:],
                               crs=self.crs, nodata_value=self.nodata)
        else:
            raise errors.NonEquivalentGridError(self, other)

    def __sub__(self, other):
        if self._equivalent_structure(other):
            return RegularGrid(copy.copy(self.transform),
                               values=self[:,:]-other[:,:],
                               crs=self.crs, nodata_value=self.nodata)
        else:
            raise errors.NonEquivalentGridError(self, other)

    def __getitem__(self, key):
        return self._bandindexer[key]

    def __setitem__(self, key, value):
        self._bandindexer[key] = value
        return

    def _equivalent_structure(self, other):
        return (self._transform == other._transform) and \
               (self.size == other.size)

    @property
    def transform(self):
        return self._transform

    @property
    def values(self):
        return self._bandindexer

    @property
    def size(self):
        return self.bands[0].size

    @property
    def nbands(self):
        return len(self.bands)

    def set_nodata_value(self, val):
        """ Redefine value used to indicate nodata.

        Parameters
        ----------
        val : number
        """
        self[:,:] = np.where(self.data_mask_full, self[:,:], val)
        self._nodata = val
        return

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

    def coordinates(self, crs=None):
        """ Return a CoordinateGenerator instance for computing the coordinates
        of indices in *crs*. This behaves like the array returned by
        `center_coords()`, but computes coordinates on demand rather than at
        once.

        Parameters
        ----------
        crs : CRS, optional
            Coordinate system of output coordinates. If None (default), output
            coordinates are not transformed from *self.crs*.

        Returns
        -------
        CoordinateGenerator

        Notes
        -----
        This is an experimental replacement for coordmesh() and center_coords().
        """
        if crs is None:
            crs = self.crs
        return CoordinateGenerator(self.transform, self.size, self.crs, crs)

    def center_coords(self):
        """ Return the coordinates of cell centers. """
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
    def origin(self):
        return self.transform[:2]

    @property
    def resolution(self):
        return self.transform[2:4]

    @property
    def skew(self):
        return self.transform[4:]

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
        ymax).

        Parameters
        ----------
        reference : string
            either 'center' or 'edge'.
        crs : karta.crs.CRS subclass
            coordinate system of output
        """
        if reference == 'center':
            x0, y0 = self.center_llref()
            n = -1
        elif reference == 'edge':
            x0, y0 = self.corner_llref()
            n = 0
        else:
            raise errors.GridError("`reference` must be 'center' or 'edge'")
        ny, nx = self.size
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
            xg0, yg0 = self.crs.transform(crs, x0, y0)
            xg1, yg1 = self.crs.transform(crs, x0, y1)
            xg2, yg2 = self.crs.transform(crs, x1, y0)
            xg3, yg3 = self.crs.transform(crs, x1, y1)
            a = min(xg0, xg1, xg2, xg3)
            b = max(xg0, xg1, xg2, xg3)
            c = min(yg0, yg1, yg2, yg3)
            d = max(yg0, yg1, yg2, yg3)
        return a, b, c, d

    def get_data_extent(self, reference='center', nodata=None, crs=None):
        """ Return the region characteristics as a tuple (xmin, xmax, ymin,
        ymax).

        Parameters
        ----------
        reference : string
            either 'center' or 'edge'.
        crs : karta.crs.CRS subclass
            coordinate system of output
        """
        if nodata is None:
            nodata = self.nodata

        if np.isnan(nodata):
            isdata = lambda a: ~np.isnan(a)
        else:
            def isdata(a):
                return a != nodata

        dx, dy = self.transform[2:4]
        sx, sy = self.transform[4:6]
        x0 = self.transform[0] + 0.5*dx + 0.5*sx
        y0 = self.transform[1] + 0.5*dy + 0.5*sy
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

        for i, row in enumerate(self[:,:]):
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
            lx, ly = self.crs.transform(crs, lx, ly)
            rx, ry = self.crs.transform(crs, rx, ry)
            bx, by = self.crs.transform(crs, bx, by)
            tx, ty = self.crs.transform(crs, tx, ty)
        return (lx, rx, by, ty)

    @property
    def data_mask_full(self):
        """ 8-bit mask of valid data cells by band """
        if np.isnan(self.nodata):
            isdata = lambda a: ~np.isnan(a)
        else:
            def isdata(a):
                return a != self.nodata
        return isdata(self[:,:])

    @property
    def data_mask(self):
        """ 8-bit mask of valid data cells, collapsed across bands """
        if len(self.bands) == 1:
            return self.data_mask_full
        else:
            return np.all(self.data_mask_full, axis=-1)

    def aschunks(self, size=(-1, -1), overlap=(0, 0), copy=True):
        """ Generator for grid chunks. This may be useful for parallel or
        memory-controlled grid processing.

        Parameters
        ----------
        size : tuple of two integers, optional
            size of the chunks to return (default approximately one-quarter of
            each dimension)
        overlap : tuple of two integers, optional
            number of pixels of overlap (default (0, 0))
        copy : bool
            whether to force returned grids to be copies
            warning: output may be a copy regardless, depending on the band class

        Yields
        ------
        RegularGrid
            smaller grids
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
            if copy:
                v = self[i0:i0+size[1], j0:j0+size[0]].copy()
            else:
                v = self[i0:i0+size[1], j0:j0+size[0]]
            yield RegularGrid(T, values=v, crs=self.crs, nodata_value=self.nodata)
            j0 += size[0]-overlap[0]

    def clip(self, xmin, xmax, ymin, ymax, crs=None):
        """ Return a clipped version of grid with cell centers constrained to a
        bounding box.

        Parameters
        ----------
        xmin : float
        xmax : float
        ymin : float
        ymax : float
        crs : karta.crs.CRS subclass, optional
        """
        if crs is not None:
            x, y = crs.transform(self.crs, [xmin, xmin, xmax, xmax],
                                           [ymin, ymax, ymin, ymax])
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

        values = self[i0:i1,j0:j1].copy()
        x0 = t[0] + j0*t[2] + i0*t[4]
        y0 = t[1] + i0*t[3] + j0*t[5]
        tnew = (x0, y0, t[2], t[3], t[4], t[5])
        return RegularGrid(tnew, values, crs=self.crs, nodata_value=self.nodata)

    def resize(self, bboxnew):
        """ Return a new grid grid with outer edges given by a bounding box.

        If the grid origin shifts by a non-integer factor of the grid
        resolution, nearest neighbour values are selected.

        Parameters
        ----------
        bboxnew : 4-tuple of floats
            (xmin, ymin, xmax, ymax).
        """
        # Define own rounding function so that half-values are treated consistently
        def _round(a):
            r = a%1
            if r <= 0.5:
                return a//1
            else:
                return a//1+1

        bb = self.bbox
        bbnew = list(bboxnew)
        dx, dy, sx, sy = self.transform[2:]

        # Redefine output box dimensions to be an integer factor of dx, dy
        bbnew[2] = bbnew[0] + dx*math.ceil((bbnew[2]-bbnew[0])/dx)
        bbnew[3] = bbnew[1] + dy*math.ceil((bbnew[3]-bbnew[1])/dy)

        ny, nx = self.size
        nxnew = int(_round((bbnew[2]-bbnew[0])/dx))
        nynew = int(_round((bbnew[3]-bbnew[1])/dy))
        Tnew = [bbnew[0], bbnew[1], dx, dy, sx, sy]

        # determine the indices of existing data on the new grid
        j0new = max(0,     int(_round((bb[0]-bbnew[0])/dx)))
        j1new = min(nxnew, int(_round((bb[2]-bbnew[0])/dx)))
        i0new = max(0,     int(_round((bb[1]-bbnew[1])/dy)))
        i1new = min(nynew, int(_round((bb[3]-bbnew[1])/dy)))

        # determine the indices of the data to copy on the old grid
        j0 = max(0,  int(_round((bbnew[0]-bb[0])/dx)))
        j1 = min(nx, int(_round((bbnew[2]-bb[0])/dx)))
        i0 = max(0,  int(_round((bbnew[1]-bb[1])/dy)))
        i1 = min(ny, int(_round((bbnew[3]-bb[1])/dy)))

        newbands = []
        for band in self.bands:
            newband = self._bndcls((nynew, nxnew), dtype=band.dtype,
                                   initval=self.nodata)
            newband[i0new:i1new, j0new:j1new] = band[i0:i1, j0:j1]
            newbands.append(newband)

        gridnew = RegularGrid(Tnew, bands=newbands, crs=self.crs,
                              nodata_value=self.nodata)
        return gridnew

    def mask_by_poly(self, polys, inplace=False):
        """ Return a grid with all elements outside the bounds of a polygon
        masked as nodata.

        Parameters
        ----------
        polys : Polygon, Multipolygon or list of Polygon instances
            region(s) defining masking boundary
        inplace : bool, optional
            if True, perform masking in-place (default False)
        """
        if getattr(polys, "_geotype", "") == "Polygon":
            polys = [polys]

        msk = None
        for poly in polys:

            x, y = poly.get_coordinate_lists(self.crs)[:2]
            if not poly.isclockwise():
                x = x[::-1]
                y = y[::-1]

            ny, nx = self.size
            _msk = mask_poly(x, y, nx, ny, self.transform)

            if msk is None:
                msk = _msk
            else:
                msk = (msk | _msk)

        if inplace:
            for band in self.bands:
                data = band[:,:]
                data[~msk] = self.nodata
                band[:,:] = data
            return self
        else:
            return RegularGrid(self.transform,
                               values=np.where(msk, self[:,:], self.nodata),
                               crs=self.crs, nodata_value=self.nodata)

    def resample(self, dx, dy, method='nearest'):
        """ Resample array to have spacing *dx*, *dy*. The grid origin remains
        in the same position.

        Parameters
        ----------
        dx : float
            cell dimension 1
        dy : float
            cell dimension 2
        method : str, optional
            interpolation method, currently 'nearest' and 'linear' supported
        """
        xmin, xmax, ymin, ymax = self.get_extent()
        if method == 'nearest':
            X, Y = np.meshgrid(np.arange(xmin, xmax, dx),
                               np.arange(ymin, ymax, dy))
            values = self.sample_nearest(X, Y)
        elif method == 'linear':
            X, Y = np.meshgrid(np.arange(xmin, xmax, dx),
                               np.arange(ymin, ymax, dy))
            values = self.sample_bilinear(X, Y)
        else:
            raise NotImplementedError('method "{0}" unavailable'.format(method))

        if values.ndim == 3:
            values = values.transpose(1, 2, 0)
        t = self._transform
        tnew = (xmin-0.5*dx-0.5*t[4], ymin-0.5*dy-0.5*t[5], dx, dy, t[4], t[5])
        return RegularGrid(tnew, values=values, crs=self.crs,
                           nodata_value=self.nodata)

    def get_positions(self, x, y):
        """ Return the float row and column indices for the point nearest
        geographical coordinates.

        Parameters
        ----------
        x, y : float or vector
            vertices of points to compute indices for

        Raises
        ------
        ValueError
            if inputs have unequal size
        """
        if hasattr(x, "__iter__"):
            x = np.array(x, dtype=np.float64)
            y = np.array(y, dtype=np.float64)
            if len(x) != len(y):
                raise ValueError("inputs must have equal size")
            return crfuncs.get_positions_vec(self.transform, x, y)
        else:
            x = np.array([x], dtype=np.float64)
            y = np.array([y], dtype=np.float64)
            I, J = crfuncs.get_positions_vec(self.transform, x, y)
            return I[0], J[0]

    def get_indices(self, x, y):
        """ Return the integer row and column indices for the point nearest
        geographical coordinates. Unlike `get_positions()`, this method raises
        and exception when points are out of range.

        Parameters
        ----------
        x, y : float or vector
            vertices of points to compute indices for

        Raises
        ------
        GridError
            points outside of Grid bbox
        """
        ny, nx = self.size
        i, j = self.get_positions(x, y)
        i = np.round(i).astype(np.int32)
        j = np.round(j).astype(np.int32)

        if hasattr(i, "__iter__"):
            if i.min() < 0 or i.max() > ny-1 or j.min() < 0 or j.max() > nx-1:
                raise errors.GridError("coordinates outside grid region "
                                       "({0})".format(self.bbox))
        else:
            if i < 0 or i > ny-1 or j < 0 or j > nx-1:
                raise errors.GridError("coordinate outside grid region "
                                       "({0})".format(self.bbox))
        return i,j

    def sample_nearest(self, x, y):
        """ Return the value nearest to coordinates using a nearest grid center
        sampling scheme.

        Parameters
        ----------
        x, y : float or vector
            vertices of points to compute indices for

        Returns
        -------
        float or vector
            The size of the output has one more dimension than the input. If
            the input are scalar, the output is a p-vector representing band
            values. If the input is a n-vector, the ouput is p x n. If the
            input is an m x n array, the output is p x m x n.

        Raises
        ------
        GridError
            points outside of Grid bbox
        """
        if not hasattr(x, "__iter__"):
            i, j = self.get_indices(x, y)
            return np.atleast_1d(self[i,j])

        if not isinstance(x, np.ndarray):
            x = np.array(x)
            y = np.array(y)

        dim = x.ndim
        I, J = self.get_indices(x.ravel(), y.ravel())

        Imn = int(np.floor(I.min())-1)
        Imx = int(np.ceil(I.max()+1))
        Jmn = int(np.floor(J.min())-1)
        Jmx = int(np.ceil(J.max()+1))

        Imn = max(Imn, 0)
        Jmn = max(Jmn, 0)
        Imx = min(Imx, self.size[0])
        Jmx = min(Jmx, self.size[1])

        I -= Imn
        J -= Jmn

        data = np.atleast_3d(self[Imn:Imx,Jmn:Jmx][I,J])
        if dim == 1:
            return data[:,:,0]
        elif dim == 2:
            m, n = x.shape
            return data.T.reshape((self.nbands, m, n))

    def sample_bilinear(self, x, y):
        """ Return the value nearest to coordinates using a bi-linear sampling
        scheme.

        Parameters
        ----------
        x, y : float or vector
            vertices of points to compute indices for

        Returns
        -------
        float or vector
            The size of the output has one more dimension than the input. If
            the input are scalar, the output is a p-vector representing band
            values. If the input is a n-vector, the ouput is p x n. If the
            input is an m x n array, the output is p x m x n.

        Raises
        ------
        GridError
            points outside of Grid bbox
        """
        if not hasattr(x, "__iter__"):
            dim = 0
            x = np.array([x])
            y = np.array([y])
        elif not isinstance(x, np.ndarray):
            x = np.array(x)
            y = np.array(y)
            dim = x.ndim
        else:
            dim = x.ndim

        I, J = self.get_positions(x.ravel(), y.ravel())

        # If the grid is large, decompressing everthing is expensive, so
        # compute only the region necessary
        if len(I) != 1:
            Imn = int(np.floor(I.min())-1)
            Imx = int(np.ceil(I.max()+1))
            Jmn = int(np.floor(J.min())-1)
            Jmx = int(np.ceil(J.max()+1))

            Imn = max(Imn, 0)
            Jmn = max(Jmn, 0)
            Imx = min(Imx, self.size[0])
            Jmx = min(Jmx, self.size[1])

            I -= Imn
            J -= Jmn
        else:
            Imn = None
            Imx = None
            Jmn = None
            Jmx = None

        data = []
        for i, band in enumerate(self.bands):
            v = self[Imn:Imx,Jmn:Jmx,i]
            if band.dtype in (np.float32, np.float64):
                data.append(crfuncs.sample_bilinear_double(I, J, v.astype(np.float64)))
            elif band.dtype in (np.int16, np.int32, np.int64):
                data.append(crfuncs.sample_bilinear_int(I, J, v.astype(np.int32)))
            elif band.dtype in (np.uint8, np.uint16, np.uint32):
                data.append(crfuncs.sample_bilinear_uint(I, J, v.astype(np.uint16)))
            else:
                raise NotImplementedError("no sample_bilinear method for dtype:"
                                          " {0}".format(band.dtype))

        if dim == 0:
            return np.array([d[0] for d in data])
        elif dim == 1:
            return np.atleast_2d(data)
        elif dim == 2:
            m, n = x.shape
            return np.concatenate([d.reshape((1, m, n)) for d in data], axis=0)

    def sample(self, *args, **kwargs):
        """ Return the values nearest positions. Positions may be:

        - a karta.Point instance
        - a karta.Multipoint instance, or
        - a pair of x, y coordinate arrays

        Parameters
        ----------
        positions : karta.Point, karta.Multipoint, or two lists of floats
            see above
        crs : karta.crs.CRS, optional
            used when coordinate lists are provided, otherwise the coordinate
            system is taken from the crs attribute of the geometry
        method : string, optional
            may be one of 'nearest', 'bilinear' (default).

        Returns
        -------
        float or vector
            The size of the output has one more dimension than the input. If
            the input are scalar, the output is a p-vector representing band
            values. If the input is a n-vector, the ouput is p x n. If the
            input is an m x n array, the output is p x m x n.

        Raises
        ------
        GridError
            points outside of Grid bbox
        """
        crs = kwargs.get("crs", None)
        method = kwargs.get("method", "bilinear")

        argerror = TypeError("`grid` takes a Point, a Multipoint, or x, y coordinate lists")
        if hasattr(args[0], "_geotype"):
            crs = args[0].crs
            if args[0]._geotype == "Point":
                x, y = args[0].get_vertex(crs=self.crs)[:2]
            elif args[0]._geotype == "Multipoint":
                x, y = args[0].get_coordinate_lists(crs=self.crs)
            else:
                raise argerror
        else:
            try:
                x = args[0]
                y = args[1]

                if crs is None:
                    crs = self.crs
                else:
                    x, y = crs.transform(self.crs, x, y)
            except IndexError:
                raise argerror
            if len(x) != len(y):
                raise argerror

        if method == "nearest":
            v = self.sample_nearest(x, y)
        elif method == "bilinear":
            v = self.sample_bilinear(x, y)
        else:
            raise ValueError("method '{0}' not available".format(method))
        return v

    def profile(self, line, resolution=None, **kw):
        """ Sample along a *Line* at a specified interval.

        Parameters
        ----------
        line : karta.Line
            defines the sampling path
        resolution : float, optional
            sample spacing, taken to be the minimum grid resolution by default

        Additional keyword arguments passed to `RegularGrid.sample` (e.g. to
        specify sampling method)

        Returns
        -------
        Multipoint
            sample points
        ndarray
            grid value at sample points
        """
        if resolution is None:
            resolution = min(self.transform[2:4])
        points = line.to_points(resolution)
        z = self.sample(*points.coordinates, **kw)
        if len(self.bands) != 1:
            for i in range(self.nbands):
                points.data.setfield("band_{0}".format(i), z[i])
        else:
            points.data.setfield("band_0", z[0])
        return points, z

    def as_warpedgrid(self):
        """ Return a copy of grid as a `WarpedGrid`. This is a more general
        grid class that has a larger memory footprint but can represent more
        flexible data layouts.
        """
        Xc, Yc = self.center_coords()
        return WarpedGrid(Xc, Yc, self[:,:].copy(), crs=self.crs)

    def to_gtiff(self, fnm, compress="PACKBITS", tiled=False, **kw):
        """ Write data to a GeoTiff file using GDAL.

        Parameters
        ----------
        fnm : str
            output file name
        compress: str or None, optional
            "PACKBITS" (default), "DEFLATE", "LZW", "LZMA", or None
        """
        return _gdal.write(fnm, self, compress=compress, **kw)

    def to_aai(self, f, reference='corner', nodata_value=-9999):
        """ Save internal data as an ASCII grid. Based on the ESRI standard,
        only isometric grids (i.e. `hdr['dx'] == hdr['dy']` can be saved.

        Parameters
        ----------
        f : str
            a file-like object or a filename
        reference : str
            specify a header reference ('center' | 'corner')
        nodata_value : number
            specify how NaNs should be represented

        Raises
        ------
        GridIOError
            input grid is not isometric and can therefore not be saved
        """
        if reference not in ('center', 'corner'):
            raise errors.GridIOError("reference in AAIGrid.tofile() must be 'center' or "
                           "'corner'")

        if np.any(self._transform[4:] != 0.0):
            raise errors.GridIOError("ESRI ASCII grids do not support skewed grids")

        ny, nx = self.bands[0].size
        x0, y0, dx, dy = self._transform[:4]
        if dx != dy:
            raise errors.GridIOError("ASCII grids require isometric grid cells")

        if not hasattr(f, 'read'):
            f = open(f, "w")

        try:
            data_a = self[:,:].copy()
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

    def to_geotiff(self, *args, **kwargs):
        return self.to_gtiff(*args, **kwargs)

    def gtiffwrite(self, *args, **kwargs):
        warnings.warn("method `gtiffwrite` has been renamed `to_gtiff`",
                FutureWarning)
        return self.to_gtiff(*args, **kwargs)

    def aaiwrite(self, *args, **kwargs):
        warnings.warn("method `aaiwrite` has been renamed `to_aai`",
                FutureWarning)
        return self.to_aai(*args, **kwargs)

class WarpedGrid(Grid):

    def __init__(self, X, Y, values, crs=None, nodata_value=None):
        """ Warped Grid class. A WarpedGrid contains a fixed number of rows and
        columns and a scalar or vector field defined on the z-axis. Grid
        spacing is not necessarily constant.

        Parameters
        ----------
        X : ndarray
            first-dimension coordinates of grid centers
        Y : ndarray
            second-dimension coordinates of grid centers
        values : ndarray
            dependent m-dimensional quantity (nrows x ncols)
        crs : karta.crs.CRS, optional
        nodata_value : number
        """

        if any(a is None for a in (X, Y, values)):
            raise errors.GridError('All of (X, Y, values) must be provided')

        if not (X.shape == Y.shape == values.shape[:2]):
            raise errors.GridError('All of (X, Y, values) must share the same '
                            'size over the first two dimensions')

        if crs is None:
            self.crs = CRS_DEFAULT
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
        """ Resample internal grid to the points defined by X, Y. """
        raise NotImplementedError

def merge(grids, weights=None):
    """ Construct a grid mosiac by averaging multiple grids. Currently limited
    to grids whose sampling is an integer translation from each other.

    Parameters
    ----------
    grids : iterable of Grid objects
        grids to combine
    weights : iterable of floats, optional
        weighting factors for computing grid averages

    Raises
    ------
    GridError, ValueError
        input grids have inconsistent transform properties
    """

    # Check grid class
    if not all(isinstance(grid, RegularGrid) for grid in grids):
        raise TypeError("all grids must be type RegularGrid")

    T = grids[0].transform
    # Check grid stretch and skew
    for i, grid in enumerate(grids[1:]):
        if grid.transform[2:6] != T[2:6]:
            raise errors.GridError("grid %d transform stretch and skew "
                    "does not match grid 1" % (i+2,))

    # Check grid offset
    excmsg = "grid %d not an integer translation from grid 1"
    for i, grid in enumerate(grids[1:]):
        if ((grid.transform[0]-T[0]) / float(T[2])) % 1 > 1e-15:
            raise ValueError(excmsg % (i+2,))
        if ((grid.transform[1]-T[1]) / float(T[3])) % 1 > 1e-15:
            raise ValueError(excmsg % (i+2,))

    # Check number of grid bands
    if len(set([len(grid.bands) for grid in grids])) != 1:
        raise ValueError("all grids to mosaic must have the same number of bands")

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
    typ = grids[0].bands[0].dtype
    outbands = []

    for iband in range(len(grids[0].bands)):
        values = np.zeros([ny, nx], dtype=typ)
        counts = np.zeros([ny, nx], dtype=np.float32)

        for grid, weight in zip(grids, normalizedweights):
            _xmin, _xmax, _ymin, _ymax = grid.get_extent(reference='edge')
            offx = int((_xmin-xmin) / T[2])
            offy = int((_ymin-ymin) / T[3])
            _ny, _nx = grid.size

            mask = grid.data_mask
            counts[offy:offy+_ny,offx:offx+_nx][mask] += weight
            _bnd = grid.bands[iband]
            values[offy:offy+_ny,offx:offx+_nx][mask] += typ(_bnd[:,:][mask]) * weight
            del mask

        validcountmask = (counts!=0.0)
        values[validcountmask] = values[validcountmask] / counts[validcountmask]
        values[~validcountmask] = grids[0].nodata

        band = CompressedBand((ny, nx), typ, initval=grids[0].nodata)
        band[:,:] = values
        outbands.append(band)

    Tmerge = [xmin, ymin] + list(T[2:])
    return RegularGrid(Tmerge, bands=outbands, crs=grids[0].crs,
                       nodata_value=grids[0].nodata)

def get_nodata(T):
    """ Return a default value for NODATA given a type

    For unsigned integer types, returns largest representable value.
    For signed integer types, returns smallest negative representable value.
    For floating point types (incl. complex), returns NaN.

    Parameters
    ----------
    T : type
        numeric type to choose NODATA for

    Raises
    ------
    ValueError
        input is not an integer, real, or complex numeric type
    """
    if T in (np.uint8, np.uint16, np.uint32, np.uint64):
        return np.iinfo(T).max
    elif issubclass(T, numbers.Integral):
        return np.iinfo(T).min
    elif issubclass(T, (numbers.Real, numbers.Complex)):
        return np.nan
    else:
        raise ValueError("No default NODATA value for type {0}".format(T))

def gridpoints(x, y, z, transform, crs):
    """ Return a grid computed by averaging point data over cells.

    Parameters
    ----------
    x : iterable
        point data x coordinates
    y : iterable
        point data y coordinates
    z : iterable
        point data values
    transform : 6-tuple of floats
        geotransform: ``[xllcorner, yllcorner, xres, yres, xskew, yskew]``
    crs : karta.crs.CRS subclass
        coordinate reference system object
    """
    ny = int((np.max(y) - transform[1]) // transform[3]) + 1
    nx = int((np.max(x) - transform[0]) // transform[2]) + 1
    grid = RegularGrid(transform,
                       values=np.zeros([ny, nx]),
                       crs=crs,
                       nodata_value=np.nan,
                       bandclass=SimpleBand)
    counts = np.zeros([ny, nx], dtype=np.int16)

    (I, J) = grid.get_indices(x, y)

    try:
        err = crfuncs.fillarray_double(grid[:,:],
                                       I.astype(np.int32),
                                       J.astype(np.int32), z, grid.nodata)
        if err != 0:
            raise RuntimeError("failure in fillarray_double")
    except ValueError:
        # Fast version works when *z* is of type double (np.float64).
        # Python fallback for other types
        for (i,j,z_) in zip(I, J, z):
            grid[i,j] += z_
            counts[i,j] += 1

        m = counts!=0
        grid[m] /= counts[m]
        grid[~m] = grid.nodata

    return grid

def mask_poly(xpoly, ypoly, nx, ny, transform):
    """ Create a grid mask based on a clockwise-oriented polygon.

    Parameters
    ----------
    xpoly, ypoly : lists of floats
        sequences of points representing polygon
    nx : int
    ny : int
        size of grid
    transform : list[float]
        affine transformation describing grid layout and origin
        ``T == [x0, y0, dx, dy, sx, sy]``
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

    # compute the grid indices of the starting point
    ta, tb, tc, td, te, tf = transform
    i0 = int(round((y0-tb - tf/tc*(x0-ta)) / (td - tf*te/tc)))
    j0 = int(round((x0-ta - te/td*(y0-tb)) / (tc - te*tf/td)))

    for el in range(1, len(xpoly)+1):
        idx = (el + i_bot) % len(xpoly)
        x1 = xpoly[idx]
        y1 = ypoly[idx]

        # (Unbounded) grid indices of the segment end points
        i1 = int(round((y1-tb - tf/tc*(x1-ta)) / (td - tf*te/tc)))
        j1 = int(round((x1-ta - te/td*(y1-tb)) / (tc - te*tf/td)))

        # If segment is horizontal or off-grid, ignore
        if ((0 <= i0 < ny) and (0 <= i1 < ny)) or (y1 != y0):

            if y1 > y0:     # upward - mark grid cells to the right

                for i in range(i0, i1):
                    if (0 <= i < ny):
                        j = int(round((i-i0) * (x1-x0)/(y1-y0) + j0))
                        if j < nx:
                            mask[i, max(0, j):] += 1

            else:           # downward - unmark grid cells to the right

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

