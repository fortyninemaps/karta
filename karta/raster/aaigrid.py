"""
Python ESRI ASCII grid driver

The primary class is AAIGrid, which can be called with a filename
argument to open a grid, or can read data from a numpy array.

Written by Nat Wilson
"""

from math import sqrt
import numpy as np
from . import grid
import traceback

class AAIGrid(grid.RegularGrid):
    """ Grid class designed to read and emit ESRI ASCII grid files. """
    _hdr = {}
    aschdr = {}
    aschdr['ncols'] = None
    aschdr['nrows'] = None
    aschdr['yllcenter'] = None
    aschdr['xllcenter'] = None
    aschdr['yllcorner'] = None
    aschdr['xllcorner'] = None
    aschdr['cellsize'] = None
    aschdr['nodata_value'] = None

    def __init__(self, incoming=None, hdr=None):
        """ Contains the table of data from an ESRI ASCII raster file.

        Parameters
        ----------
        incoming : either a filename or an ndarray containing the
                   raster data. If incoming is not provided, then
                   read() must be called later with a suitable
                   filename.
        hdr : dictionary containing AAIGrid header information

        Header dictionaries contain the following keys:
            ncols           (int)
            nrows           (int)
            either:
                xllcenter   (int)
                yllcenter   (int)
            or:
                xllcorner   (int)
                yllcorner   (int)
            cellsize        (int)
            nodata_value    (int)
        """

        if hdr is not None:
            self.aschdr = hdr

        if incoming is not None:
            try:
                # assume incoming is a file
                self.fromfile(incoming)
            except (TypeError, UnicodeDecodeError):
                # isn't a file, so hopefully an ndarray
                self.fromarray(incoming, hdr)
        elif hdr is not None:
            self.data = np.empty((hdr['nrows'], hdr['ncols']))
        else:
            self.data = None

        self._check_header()
        self._enforce_hdr_consistency()
        return

    def __str__(self):
        if self.data is None:
            stat = 'empty'
        else:
            stat = 'populated'
        return "{status} AAIGrid instance".format(status=stat)

    def __len__(self):
        if self.data is None:
            return 0
        else:
            return len(self.data)

    def __add__(self, other):
        if isinstance(other, AAIGrid):
            if (self.data is not None) and (other.data is not None):
                return AAIGrid(self.data + other.data, self.aschdr)
            else:
                raise ValueError("cannot add NoneType data array")
        else:
            try:
                return AAIGrid(self.data + other, self.aschdr)
            except AttributeError:
                raise AAIError("self.data not defined")
            except Exception:
                raise AAIError("addition with type {0} not "
                               "defined".format(type(other)))

    def __sub__(self, other):
        if isinstance(other, AAIGrid):
            if (self.data is not None) and (other.data is not None):
                return AAIGrid(self.data - other.data, self.aschdr)
            else:
                raise ValueError("cannot subtract NoneType data array")
        else:
            try:
                return AAIGrid(self.data - other, self.aschdr)
            except AttributeError:
                raise AAIError("self.data not defined")
            except Exception:
                raise AAIError("subtraction with type {0} not "
                               "defined".format(type(other)))

    def __mul__(self, other):
        if isinstance(other, AAIGrid):
            if (self.data is not None) and (other.data is not None):
                return AAIGrid(self.data * other.data, self.aschdr)
            else:
                raise ValueError("cannot multiply NoneType data array")
        else:
            try:
                return AAIGrid(self.data * other, self.aschdr)
            except AttributeError:
                raise AAIError("self.data not defined")
            except:
                raise AAIError("multiplication with type{0} not "
                                "defined".format(type(other)))

    def __div__(self, other):
        if isinstance(other, AAIGrid):
            if (self.data is not None) and (other.data is not None):
                return AAIGrid(self.data / other.data, self.aschdr)
            else:
                raise ValueError("cannot divide NoneType data array")
        else:
            try:
                return AAIGrid(self.data / other, self.aschdr)
            except AttributeError:
                raise AAIError("self.data not defined")
            except:
                raise AAIError("division with type{0} not "
                                "defined".format(type(other)))

    def _check_header(self, hdr=None):
        """ Make sure that all required header records are present. """

        if hdr is None:
            hdr = self.aschdr

        for field in ('ncols', 'nrows', 'cellsize'):
            if hdr.get(field) is None:
                raise AAIError("{0} not defined".format(field.upper()))

        for field in ('ncols', 'nrows'):
            hdr[field] = int(hdr[field])

        if hdr.get('nodata_value') is None:
            hdr['nodata_value'] = -9999

        hdr = self._check_header_references(hdr)

        return hdr

    def _check_header_references(self, hdr):
        """ Make sure that both corner and center reference coordinates are
        present. """
        d = hdr['cellsize']
        if hdr.get('yllcenter') is None:
            if hdr.get('yllcorner') is None:
                raise AAIError('y reference not defined')
            else:
                hdr['yllcenter'] = hdr['yllcorner'] + d / 2.0
        else:
            hdr['yllcorner'] = hdr['yllcenter'] - d / 2.0

        if hdr.get('xllcenter') is None:
            if hdr.get('xllcorner') is None:
                raise AAIError('x reference not defined')
            else:
                hdr['xllcenter'] = hdr['xllcorner'] + d / 2.0
        else:
            hdr['xllcorner'] = hdr['xllcenter'] - d / 2.0
        return hdr

    def _enforce_hdr_consistency(self):
        """ Make sure that base `hdr` remains consistent with AAIGrid-specific
        `aschdr`. """
        self._hdr['dx'] = self.aschdr['cellsize']
        self._hdr['dy'] = self.aschdr['cellsize']
        self._hdr['nx'] = self.aschdr['ncols']
        self._hdr['ny'] = self.aschdr['nrows']
        self._hdr['xllcorner'] = self.aschdr['xllcorner']
        self._hdr['yllcorner'] = self.aschdr['yllcorner']
        return

    def get_indices(self, x, y):
        """ Return the column and row indices for the point nearest
        geographical coordinates (x, y). """
        if self.data is None:
            raise AAIError('no raster to query')

        self.aschdr = self._check_header_references(self.aschdr)

        x0 = self.aschdr['xllcenter']
        y0 = self.aschdr['yllcenter']
        d = self.aschdr['cellsize']
        nx, ny = self.data.shape

        xi = np.clip(
                (np.around((np.array(x) - x0) / d)).astype(int),
                0, nx-1)
        yi = np.clip(
                self.data.shape[0] -
                    (np.around((np.array(y)-y0) / d)).astype(int) - 1,
                0, ny-1)
        return xi, yi

    def get_region(self, reference='center'):
        """ Return the region characteristics as a tuple, relative to either
        cell centers or the grid edges. """
        if reference == 'center':
            x0, y0  = self.aschdr['xllcenter'], self.aschdr['yllcenter']
            n = 0
        elif reference == 'edge':
            x0, y0  = self.aschdr['xllcenter'], self.aschdr['yllcenter']
            n = 1
        else:
            raise AAIError("`reference` must be 'center' or 'edge'")
        return (x0, x0 + (self.aschdr['ncols'] + n) * self.aschdr['cellsize'],
                y0, y0 + (self.aschdr['nrows'] + n) * self.aschdr['cellsize'])

    def coordmesh(self, grid_registration='center'):
        """ Return a pair of arrays containing the *X* and *Y* coordinates of
        the grid. """
        if grid_registration == 'center':
            xll = self.aschdr['xllcenter']
            yll = self.aschdr['yllcenter']
        elif grid_registration == 'corner':
            xll = self.aschdr['xllcorner']
            yll = self.aschdr['yllcorner']
        else:
            raise AAIError("grid_registration must be 'center' or 'corner'\n")

        X = (xll + np.arange(self.data.shape[1]) * self.aschdr['cellsize'])
        Y = (yll + np.arange(self.data.shape[0])[::-1] * self.aschdr['cellsize'])
        return np.meshgrid(X, Y)

    def max(self):
        """ Return the maximum non-nan in self.data. """
        return self.data[np.isnan(self.data)==False].max()

    def min(self):
        """ Return the minimum non-nan in self.data. """
        return self.data[np.isnan(self.data)==False].min()

    def read(self, fnm):
        """ Read an existing ASCII grid file. """
        return self.fromfile(fnm)

    def fromfile(self, fnm):
        """ Read an existing ASCII grid file. """
        try:
            h = []
            with open(fnm, 'r') as f:
                # Count the number of header entries
                cnt = 0
                while True:
                    l = f.readline()
                    if l.split(None, 1)[0].lower() in ['nrows', 'ncols',
                            'yllcenter', 'xllcenter', 'yllcorner', 'xllcorner',
                            'cellsize', 'nodata_value']:
                        cnt += 1
                    else:
                        break

                # Read the header, then the array data
                f.seek(0)
                for _ in range(cnt):
                    h.append(f.readline())
                data = f.readlines()

        except IOError:
            raise AAIError('error while trying to open {0}'.format(fnm))

        # Store data internally
        try:
            hs = [rec.split(None, 1) for rec in h]
            hdr = dict([(rec[0].lower(), float(rec[1])) for rec in hs])
            if 'yllcenter' not in hdr.keys():
                hdr['yllcenter'] = None

            if 'xllcenter' not in hdr.keys():
                hdr['xllcenter'] = None

            if 'yllcorner' not in hdr.keys():
                hdr['yllcorner'] = None

            if 'xllcorner' not in hdr.keys():
                hdr['xllcorner'] = None

            if 'nodata_value' not in hdr.keys():
                hdr['nodata_value'] = -9999

        except:
            traceback.print_exc()
            return

        try:
            hdr = self._check_header(hdr)
        except AAIError as s:
            print(s)
            return

        self.aschdr = hdr

        try:
            f = lambda l: [float(i) for i in l.split()]
            data_f = map(f, data)
            data_a = np.array(data_f)
            data_a[data_a==hdr['nodata_value']] = np.nan
            self.data = data_a
        except:
            traceback.print_exc()

        return

    def fromarray(self, A, hdr):
        """ Read data from a 2d numpy array. The hdr argument is a dictionary
        of header information.

        ncols : int
            Number of columns

        nrows : int
            Number of rows

        cellsize : float
            Cell dimension

        nodata_value : float
            Value for missing data

        xllcenter OR xllcorner : float
            Lower left cell *x*-registration

        yllcenter OR yllcorner : float
            Lower left cell *y*-registration
        """
        # Should add sanity checks
        self.aschdr = hdr
        #A[np.isnan(A)] = hdr['nodata_value']
        self.data = A.copy()[:,:]

    def tofile(self, f, reference='center'):
        """ Save internal data to *f*, either a file-like object or a filename.

        Parameters
        ----------
        f : filename or file handle
        reference : grid registration ("center" or "corner")
        """
        if self.data is None:
            raise AAIError("no data to write!")

        if reference not in ('center', 'corner'):
            raise AAIError("reference in AAIGrid.tofile() must be 'center' or "
                           "'corner'")

        try:
            f.read()
        except AttributeError:
            f = open(f, "w")
        except IOError:
            pass

        try:
            data_a = self.data.copy()
            data_a[np.isnan(data_a)] = self.aschdr['nodata_value']

            f.write("NCOLS {0}\n".format(self.aschdr['ncols']))
            f.write("NROWS {0}\n".format(self.aschdr['nrows']))
            if reference == 'center':
                f.write("XLLCENTER {0}\n".format(self.aschdr['xllcenter']))
                f.write("YLLCENTER {0}\n".format(self.aschdr['yllcenter']))
            elif reference == 'corner':
                f.write("XLLCORNER {0}\n".format(self.aschdr['xllcorner']))
                f.write("YLLCORNER {0}\n".format(self.aschdr['yllcorner']))
            f.write("CELLSIZE {0}\n".format(self.aschdr['cellsize']))
            f.write("NODATA_VALUE {0}\n".format(self.aschdr['nodata_value']))
            f.writelines([str(row).replace(',','')[1:-1] +
                            '\n' for row in data_a.tolist()])

        except:
            traceback.print_exc()
        finally:
            f.close()
        return

    def toarray(self):
        """ Return a 2d array with internal data. Header information can be
        obtained with AAIGrid.hdr.
        """
        return self.data.copy()[::-1,:]

    def minmax(self):
        """ Return the minimum and maximum value of data array. """
        A = self.data[np.isnan(self.data)==False]
        return (A.min(), A.max())

    def sample(self, x, y):
        """ Return the value nearest to *x*, *y*. """
        xi, yi = self.get_indices(x, y)
        ncols = self.aschdr['ncols']
        nrows = self.aschdr['nrows']
        if (np.any(xi < 0) or np.any(xi >= ncols) or
            np.any(yi < 0) or np.any(yi >= nrows)):
            raise AAIError("coordinates are outside grid region",
                detail="({0}, {1}), ({2}, {3})".format(xi, yi, ncols, nrows))
        else:
            z = self.data[yi, xi]
        return z

    def get_profile(self, segments, resolution=10.0):
        """ Sample along a line defined as `segments`. Does not interpolate.

        Parameters
        ----------
        segments : iterable containing (x,y) pairs
        resolution : sample spacing

        Returns
        -------
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
                z.append(self.data[yi, xi])
                p += resolution
            p -= d

        return np.array(z)

    def clip(self, bounds):
        """ Clip the z-range in place to bounds = [min, max]. """
        self.data = self.data.clip(bounds[0], bounds[1])
        return

    def resample(self, cellsize, method='nearest'):
        """ Resample array to have `cellsize`. Modify header in-place.

        Parameters
        ----------
        cellsize : cell dimension, float
        method : interpolation method, string ('nearest', 'linear')
        """
        xllcenter = self.aschdr['xllcenter']
        yllcenter = self.aschdr['yllcenter']
        xurcenter = xllcenter + self.aschdr['cellsize'] * self.aschdr['ncols']
        yurcenter = yllcenter + self.aschdr['cellsize'] * self.aschdr['nrows']
        nx = int((xurcenter - xllcenter) // cellsize)
        ny = int((yurcenter - yllcenter) // cellsize)
        dimratio = cellsize / self.aschdr['cellsize']

        if method == 'nearest':
            JJ, II = np.meshgrid(np.arange(nx), np.arange(ny))
            srcII = np.around(II * dimratio) \
                            .astype(int) \
                            .clip(0, self.aschdr['nrows'] - 1)
            srcJJ = np.around(JJ * dimratio) \
                            .astype(int).\
                            clip(0, self.aschdr['ncols'] - 1)
            self.data = self.data[srcII, srcJJ]
        else:
            raise NotImplementedError('method "{0}" not '
                                      'implemented'.format(method))

        self.aschdr['cellsize'] = float(cellsize)
        self.aschdr['nrows'] = ny
        self.aschdr['ncols'] = nx
        self.aschdr['xllcorner'] = self.aschdr['xllcenter'] - (0.5 * cellsize)
        self.aschdr['yllcorner'] = self.aschdr['yllcenter'] - (0.5 * cellsize)
        self._enforce_hdr_consistency()
        return

    def resize(self, te):
        """ Resize array to fit within extent given by te. If the new
        dimensions are smaller, the data is clipped. If they are larger,
        nan padding is added.

        Parameters
        ----------
        te : tuple of center coordinates in the form (xmin, xmax, ymin, ymax).

        Returns None.
        """
        if self.data is not None:
            xmin1 = self.aschdr['xllcenter']
            xmax1 = xmin1 + self.aschdr['cellsize'] * (self.aschdr['ncols'] - 1)
            ymin1 = self.aschdr['yllcenter']
            ymax1 = ymin1 + self.aschdr['cellsize'] * (self.aschdr['nrows'] - 1)

            xmin2 = te[0]
            xmax2 = te[1]
            ymin2 = te[2]
            ymax2 = te[3]

            data_a = self.data.copy()
            ny, nx = data_a.shape

            # The left side
            dx = int(np.floor((xmin2-xmin1) / float(self.aschdr['cellsize'])))
            if dx < 0:
                data_a = np.hstack([np.full((ny, dx), np.nan), data_a])
            elif dx > 0:
                data_a = data_a[:,dx:]
            self.aschdr['xllcenter'] = xmin1 + self.aschdr['cellsize'] * dx

            # The right side
            dx = int(np.floor((xmax2-xmax1) / float(self.aschdr['cellsize'])))
            if dx > 0:
                data_a = np.hstack([data_a, np.full((ny, dx), np.nan)])
            elif dx < 0:
                data_a = data_a[:,:dx]
            self.aschdr['ncols'] = data_a.shape[1]

            ny, nx = data_a.shape

            # The bottom
            dy = int(np.ceil((ymin2-ymin1) / float(self.aschdr['cellsize'])))
            if dy < 0:
                data_a = np.vstack([data_a, np.full((dy, nx), np.nan)])
            elif dy > 0:
                data_a = data_a[dy:,:]
            self.aschdr['yllcenter'] = ymin1 + self.aschdr['cellsize'] * dy

            # The top
            dy = int(np.floor((ymax2-ymax1) / float(self.aschdr['cellsize'])))
            if dy > 0:
                data_a = np.vstack([np.full((dy, nx), np.nan), data_a])
            elif dy < 0:
                data_a = data_a[:dy,:]
            self.aschdr['nrows'] = data_a.shape[0]

            self.data = data_a
            self._check_header()
            self._enforce_hdr_consistency()

        else:
            raise AAIError("no data to resize")
        return

class AAIError(Exception):
    """ Exceptions related to ESRI ASCII grid driver. """
    def __init__(self, value, detail=None):
        self.value = value
        self.detail = detail
    def __str__(self):
        if self.detail is not None:
            s = self.value + ': ' + self.detail
        else:
            s = self.value
        return repr(s)

