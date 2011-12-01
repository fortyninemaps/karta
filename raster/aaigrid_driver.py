"""
Python ESRI ASCII grid driver

The primary class is AAIGrid, which can be called with a filename
argument to open a grid, or can read data from a numpy array.

Written by Nat Wilson
"""

from math import sin, cos, sqrt
import numpy as np
import traceback

class AAIGrid(object):
    def __init__(self, incoming=None, hdr=None):
        """ Contains the table of data from an ESRI ASCII raster file.

        *incoming*  : either a filename or an ndarray containing the
                      raster data. If incoming is not provided, then
                      read() must be called later with a suitable
                      filename.

        *hdr*       : dictionary containing AAIGrid header information
        """

        self.data = None
        self.hdr = {}
        self.hdr['ncols'] = None
        self.hdr['nrows'] = None
        self.hdr['yllcenter'] = None
        self.hdr['xllcenter'] = None
        self.hdr['yllcorner'] = None
        self.hdr['xllcorner'] = None
        self.hdr['cellsize'] = None
        self.hdr['nodata_value'] = None
        if incoming is not None:
            try:
                # assume incoming is a file
                self.fromfile(incoming)
            except TypeError:
                # isn't a file, so hopefully an ndarray
                self.fromarray(incoming, hdr)

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
                return AAIGrid(self.data + other.data, self.hdr)
            else:
                raise ValueError("cannot add NoneType data array")
        else:
            raise AAIError("addition with other types not defined")

    def __sub__(self, other):
        if isinstance(other, AAIGrid):
            if (self.data is not None) and (other.data is not None):
                return AAIGrid(self.data - other.data, self.hdr)
            else:
                raise ValueError("cannot subtract NoneType data array")
        else:
            raise AAIError("subtraction with other types not defined")

    def __mul__(self, other):
        if isinstance(other, AAIGrid):
            if (self.data is not None) and (other.data is not None):
                return AAIGrid(self.data * other.data, self.hdr)
            else:
                raise ValueError("cannot multiply NoneType data array")
        else:
            raise AAIError("multiplication with other types not defined")

    def __div__(self, other):
        if isinstance(other, AAIGrid):
            if (self.data is not None) and (other.data is not None):
                return AAIGrid(self.data / other.data, self.hdr)
            else:
                raise ValueError("cannot divide NoneType data array")
        else:
            raise AAIError("division with other types not defined")

    def _check_header(self, hdr):
        """ Make sure that all required header records are present. """
        if hdr['ncols'] is None:
            raise AAIError('NCOLS not defined')
        else:
            hdr['ncols'] = int(hdr['ncols'])

        if hdr['nrows'] is None:
            raise AAIError('NROWS not defined')
        else:
            hdr['nrows'] = int(hdr['nrows'])

        if hdr['cellsize'] is None:
            raise AAIError('CELLSIZE not defined')

        if hdr['nodata_value'] is None:
            hdr['nodata_value'] = -9999

        if (hdr['yllcenter'] is None):
            if (hdr['yllcorner'] is None):
                raise AAIError('y reference not defined')
            else:
                hdr['yllcenter'] = hdr['yllcorner'] + hdr['cellsize'] / 2.0
        else:
            hdr['yllcorner'] = hdr['yllcenter'] - hdr['cellsize'] / 2.0

        if (hdr['xllcenter'] is None):
            if (hdr['xllcorner'] is None):
                raise AAIError('x reference not defined')
            else:
                hdr['xllcenter'] = hdr['xllcorner'] + hdr['cellsize'] / 2.0
        else:
            hdr['xllcorner'] = hdr['xllcenter'] - hdr['cellsize'] / 2.0
        return hdr

    def _get_indices(self, x, y):
        """ Return the column and row indices for the point nearest
        geographical coordinates (x, y). """
        if self.data is None:
            raise AAIError('no raster to query')
            return None
        d = self.hdr['cellsize']
        if self.hdr['xllcenter'] is None:
            if self.hdr['xllcorner'] is not None:
                self.hdr['xllcenter'] = self.hdr['xllcorner'] + d/2.0
            else:
                raise AAIError('x reference not defined')
        if self.hdr['yllcenter'] is None:
            if self.hdr['yllcorner'] is not None:
                self.hdr['yllcenter'] = self.hdr['yllcorner'] + d/2.0
            else:
                raise AAIError('y reference not defined')
        x0 = self.hdr['xllcenter']
        y0 = self.hdr['yllcenter']
        xmax = x0 + d * self.hdr['ncols']
        ymax = y0 + d * self.hdr['nrows']
        xi = int(round((x-x0) / d))
        yi = self.data.shape[0] - int(round((y-y0) / d)) - 1
        return xi, yi

    def max(self):
        """ Return the maximum non-nan in self.data. """
        return self.data[np.isnan(self.data)==False].max()

    def min(self):
        """ Return the minimum non-nan in self.data. """
        return self.data[np.isnan(self.data)==False].min()

    def Read(self, fnm):
        """ Read an existing ASCII grid file. """
        return self.fromfile(fnm)

    def read(self, fnm):
        """ Read an existing ASCII grid file. """
        return self.fromfile(fnm)

    def FromFile(self, fnm):
        """ Read an existing ASCII grid file. """
        return self.fromfile(fnm)

    def fromfile(self, fnm):
        """ Read an existing ASCII grid file. """
        try:
            h = []
            # Count the number of header entries
            with open(fnm, 'r') as f:
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
            with open(fnm, 'r') as f:
                for i in range(cnt):
                    h.append(f.readline())
                data = f.readlines()
        except IOError:
            raise AAIError('error while trying to open {0}'.format(fnm))
            return
        # Store data internally
        try:
            hs = [rec.split(None,1) for rec in h]
            hdr = dict([(rec[0].lower(), float(rec[1])) for rec in hs])
            if 'yllcenter' not in hdr.keys():
                hdr['yllcenter'] = None

            if 'xllcenter' not in hdr.keys():
                hdr['xllcenter'] = None

            if 'yllcorner' not in hdr.keys():
                hdr['yllcorner'] = None

            if 'xllcorner' not in hdr.keys():
                hdr['xllcorner'] = None

        except:
            traceback.print_exc()
            return

        try:
            hdr = self._check_header(hdr)
        except AAIError as s:
            print s
            return

        self.hdr = hdr

        try:
            f = lambda l: [float(i) for i in l.split()]
            data_f = map(f, data)
            data_a = np.array(data_f)
            data_a[data_a==hdr['nodata_value']] = np.nan
            self.data = data_a
        except:
            traceback.print_exc()

        return

    def FromArray(self, A, hdr):
        return self.fromarray(A, hdr)

    def fromarray(self, A, hdr):
        """ Read data from a 2d numpy array. The hdr argument is a dictionary
        of header information, including:
            ncols: [int]
            nrows: [int]
            cellsize: [float]
            nodata_value: [float]
            xllcenter OR xllcorner: [float]
            yllcenter OR yllcorner: [float]
        """
        # Should add sanity checks
        self.hdr = hdr
        #A[np.isnan(A)] = hdr['nodata_value']
        self.data = A.copy()[:,:]
        pass

    def ToFile(self, f, reference='center'):
        return self.tofile(f, reference='center')

    def tofile(self, f, reference='center'):
        """ Save internal data to f, either a file-like object or a
        filename. """
        if self.data is None:
            raise AAIError("no data to write!")

        if reference not in ('center', 'corner'):
            raise AAIError("reference in AAIGrid.tofile() must be 'center' or 'corner'")

        try:
            f.read()
        except AttributeError:
            f = open(f, "w")
        except IOError:
            pass

        try:
            data_a = self.data.copy()
            data_a[np.isnan(data_a)] = self.hdr['nodata_value']

            f.write("NCOLS {0}\n".format(self.hdr['ncols']))
            f.write("NROWS {0}\n".format(self.hdr['nrows']))
            if reference == 'center':
                f.write("XLLCENTER {0}\n".format(self.hdr['xllcenter']))
                f.write("YLLCENTER {0}\n".format(self.hdr['yllcenter']))
            elif reference == 'corner':
                f.write("XLLCORNER {0}\n".format(self.hdr['xllcorner']))
                f.write("YLLCORNER {0}\n".format(self.hdr['yllcorner']))
            f.write("CELLSIZE {0}\n".format(self.hdr['cellsize']))
            f.write("NODATA_VALUE {0}\n".format(self.hdr['nodata_value']))
            f.writelines([str(row).replace(',','')[1:-1] +
                            '\n' for row in data_a.tolist()])

        except:
            traceback.print_exc()
        finally:
            f.close()
        return

    def ToArray(self):
        return self.toarray()

    def toarray(self):
        """ Return a 2d array with internal data. Header information can be
        obtained with AAIGrid.hdr.
        This is just a nicer way of calling AAIGrid.data.copy().
        """
        return self.data.copy()[::-1,:]

    def MinMax(self):
        return self.minmax()

    def minmax(self):
        """ Return the minimum and maximum value of data array. """
        A = self.data[np.isnan(self.data)==False]
        return (A.min(), A.max())

    def Sample(self, x, y):
        return self.sample(x, y)

    def sample(self, x, y):
        """ Return the value nearest to x, y, as well as center coordinates of
        the grid cell actually sampled. """
        xi, yi = self._get_indices(x, y)
        d = self.hdr['cellsize']
        ncols = self.hdr['ncols']
        nrows = self.hdr['nrows']
        if (xi < 0) or (xi >= ncols) or (yi < 0) or (yi >= nrows):
            raise AAIError("coordinates are outside grid region",
                detail="({0}, {1}), ({2}, {3})".format(xi, yi, ncols, nrows))
            return
        else:
            z = self.data[yi, xi]
            ys = yi * self.hdr['cellsize'] + self.hdr['yllcenter']
            xs = xi * self.hdr['cellsize'] + self.hdr['xllcenter']
            return z, (xs, ys)

    def GetProfile(self, x0, y0, xf, yf, resolution=10.0):
        return self.get_profile(x0, y0, xf, yf, resolution=10.0)

    def get_profile(self, x0, y0, xf, yf, resolution=10.0):
        """ Sample along a line defined as *segments*.

        *segments*      :   iterable containing (x,y) pairs
        *resolution*    :   sample spacing

        Returns an ndarray.
        Does not interpolate.
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
                xi, yi = self._get_indices(x0+fx, y0+fy)
                z.append(self.data[yi, xi])
                p += resolution
            p -= d

        return np.array(z)

    def Clip(self, bounds):
        return self.clip(bounds)

    def clip(self, bounds):
        """ Clip the z-range in place to bounds = [min, max]. """
        self.data = self.data.clip(bounds[0], bounds[1])
        return

    def Resize(self, te):
        return self.resize(te)

    def resize(self, te):
        """ Resize array to fit within extents given by te. If the new
        dimensions are smaller, the data is clipped. If they are larger,
        nan padding is added.

        *te*        : tuple in the form (xmin, xmamx, ymin, ymax).

        Returns None.
        """
        if self.data is not None:
            xmin1 = self.hdr['xllcenter']
            xmax1 = xmin1 + self.hdr['cellsize'] * (self.hdr['ncols'] - 1)
            ymin1 = self.hdr['yllcenter']
            ymax1 = ymin1 + self.hdr['cellsize'] * (self.hdr['nrows'] - 1)

            xmin2 = te[0]
            xmax2 = te[1]
            ymin2 = te[2]
            ymax2 = te[3]

            data_a = self.data.copy()
            ny, nx = data_a.shape

            if xmin2 < xmin1:
                dx = int(np.ceil((xmin2-xmin1) / float(self.hdr['cellsize'])))
                print dx
                data_a = np.hstack([np.nan*np.ones([ny, -dx]), data_a])
            elif xmin2 > xmin1:
                dx = int(np.ceil((xmin2-xmin1) / float(self.hdr['cellsize'])))
                data_a = data_a[:,dx:]
            else:
                dx = 0
            self.hdr['xllcenter'] = xmin1 + self.hdr['cellsize'] * dx

            if xmax2 > xmax1:
                dx = int(np.floor((xmax2-xmax1) / float(self.hdr['cellsize'])))
                print dx
                data_a = np.hstack([data_a, np.nan*np.ones([ny, dx])])
            elif xmax2 < xmax1:
                dx = int(np.floor((xmin2-xmin1) / float(self.hdr['cellsize'])))
                data_a = data_a[:,:dx]
            else:
                dx = 0
            self.hdr['ncols'] = data_a.shape[1]

            ny, nx = data_a.shape

            if ymin2 < ymin1:
                dy = int(np.ceil((ymin2-ymin1) / float(self.hdr['cellsize'])))
                print dy
                data_a = np.vstack([data_a, np.nan*np.ones([-dy, nx])])
            elif ymin2 > ymin1:
                dy = int(np.ceil((ymin2-ymin1) / float(self.hdr['cellsize'])))
                data_a = data_a[dy:,:]
            else:
                d = 0
            self.hdr['yllcenter'] = ymin1 + self.hdr['cellsize'] * dy

            if ymax2 > ymax1:
                dy = int(np.floor((ymax2-ymax1) / float(self.hdr['cellsize'])))
                print dy
                data_a = np.vstack([np.nan*np.ones([dy, nx]), data_a])
            elif ymax2 < ymax1:
                dy = int(np.floor((ymin2-ymin1) / float(self.hdr['cellsize'])))
                data_a = data_a[:dy,:]
            else:
                dy = 0
            self.hdr['nrows'] = data_a.shape[0]

            self.data = data_a
            self.hdr = self._check_header(self.hdr)

        else:
            raise AAIError("no data to resize")
        return

Grid = AAIGrid      # For backwards compatibility

class AAIError(Exception):
    def __init__(self, value, detail=None):
        self.value = value
        self.detail = detail
    def __str__(self):
        if self.detail is not None:
            s = self.value + ': ' + self.detail
        else:
            s = self.value
        return repr(s)
