""" IO interface to GeoTiffs using GDAL. """

import sys
import numpy as np
from .. import errors

try:
    import osgeo.gdal
    import osgeo.osr
    import osgeo.gdalconst as gc
    osgeo.gdal.UseExceptions()
    HASGDAL = True
except ImportError:
    HASGDAL = False

class GdalBandArrayInterface(object):
    """ Imitates an ndarray well-enough to back a Grid instance, but reads data
    from an disk-bound datasource """

    def __init__(self, band, dataset):
        self.band = band
        self.dataset = dataset
        return

    def __del__(self):
        self.dataset = None
        self.band = None

    def __getitem__(self, idx):
        if isinstance(idx, tuple):
            if isinstance(idx[0], slice):
                yslc, xslc = idx
            else:
                yslc = slice(idx[0], idx[0]+1, 1)
                xslc = slice(idx[1], idx[1]+1, 1)
        else:
            yslc = idx
            xslc = slice(0, None)

        ny, nx = self.shape
        xstart, xend, xstep = xslc.indices(nx)
        ystart, yend, ystep = yslc.indices(ny)
        x0 = min(xstart, xend)
        y0 = min(ny-ystart, ny-yend)
        ret = self.band.ReadAsArray(x0, y0,
                                    abs(xend-xstart), abs(yend-ystart))
        if xstep < 0:
            ret = ret[:,::-1]

        if ystep > 0:       # Flip from GDAL convention to karta convention
            ret = ret[::-1]

        if xstep != 1:
            ret = ret[::abs(ystep)]

        if ystep != 1:
            ret = ret[:,::abs(xstep)]
        return ret

    @property
    def shape(self):
        return (self.dataset.RasterYSize, self.dataset.RasterXSize)

    @property
    def dtype(self):
        return np.dtype(numpy_dtype(self.band.DataType))


def SRS_from_WKT(s):
    """ Return Proj.4 string, semimajor axis, and flattening """
    sr = osgeo.osr.SpatialReference()
    sr.ImportFromWkt(s)
    return sr

def numpy_dtype(dt_int):
    """ Return a numpy dtype that matches the band data type. """
    name = osgeo.gdal.GetDataTypeName(dt_int)
    if name == "Byte":
        return np.uint8
    elif name == "UInt16":
        return np.uint16
    elif name == "Int16":
        return np.int16
    elif name == "UInt32":
        return np.uint32
    elif name == "Int32":
        return np.int32
    elif name == "Float32":
        return np.float32
    elif name == "Float64":
        return np.float64
    elif name in ("CInt16", "CInt32", "CFloat32", "CFloat64"):
        return np.complex64
    else:
        raise TypeError("GDAL data type {0} unknown to karta".format(dt_int))

def gdal_type(dtype):
    """ Return a GDAL type that matches numpy dtype """
    if dtype == np.uint16:
        return osgeo.gdal.GDT_UInt16
    elif dtype == np.int16:
        return osgeo.gdal.GDT_Int16
    elif dtype == np.int32:
        return osgeo.gdal.GDT_Int32
    elif dtype == np.int32:
        return osgeo.gdal.GDT_Int32
    elif dtype == np.float32:
        return osgeo.gdal.GDT_Float32
    elif dtype == np.float64:
        return osgeo.gdal.GDT_Float64
    elif dtype == np.complex64:
        return osgeo.gdal.GDT_CFloat64
    else:
        raise TypeError("GDAL equivalent to type {0} unknown".format(dtype))

def read(fnm, band, in_memory):
    """ Read a GeoTiff file and return a numpy array and a dictionary of header
    information.
    
    fnm: input datasource
    band: band number (1...)
    in_memory: boolean indicating whether array should be read fully into memory

    Returns and array-like object and a dictionary of metadata
    """
    if not HASGDAL:
        raise errors.MissingDependencyError("requires osgeo.gdal")

    hdr = dict()
    dataset = osgeo.gdal.Open(fnm, gc.GA_ReadOnly)

    try:
        hdr["nx"] = dataset.RasterXSize
        hdr["ny"] = dataset.RasterYSize

        transform = dataset.GetGeoTransform()
        if transform is not None:
            hdr["dx"] = transform[1]
            hdr["dy"] = transform[5]
            hdr["xulcorner"] = transform[0]
            hdr["yulcorner"] = transform[3]
            hdr["sx"] = transform[2]
            hdr["sy"] = transform[4]
        else:
            raise AttributeError("No GeoTransform in geotiff file")

        sr = SRS_from_WKT(dataset.GetProjectionRef())
        hdr["srs"] = {"proj4": sr.ExportToProj4(),
                      "semimajor": sr.GetSemiMajor(),
                      "flattening": sr.GetInvFlattening()}

        max_dtype = 0
        rasterband = dataset.GetRasterBand(band)
        nx = rasterband.XSize
        ny = rasterband.YSize
        if rasterband.DataType > max_dtype:
            max_dtype = rasterband.DataType

        if in_memory:
            arr = np.empty((hdr["ny"], hdr["nx"]), dtype=numpy_dtype(max_dtype))
            dtype = numpy_dtype(rasterband.DataType)
            arr[:,:] = rasterband.ReadAsArray(buf_obj=np.empty([ny, nx], dtype=dtype))
        else:
            arr = GdalBandArrayInterface(rasterband, dataset)

    finally:
        if in_memory:
            rasterband = None
            dataset = None
    return arr, hdr

def srs_from_crs(crs):
    srs = osgeo.osr.SpatialReference()
    # SpatialReference can't parse 'lonlat'
    proj4 = crs.get_proj4().replace("lonlat", "latlong")
    srs.ImportFromProj4(proj4)
    return srs

def write(fnm, grid, compress=None, **kw):
    """ Write a grid-like object to *fnm* """
    if not HASGDAL:
        raise errors.MissingDependencyError("requires osgeo.gdal")

    co = []
    if compress == "LZW":
        co.append("COMPRESS=LZW")
    elif compress == "PACKBITS":
        co.append("COMPRESS=PACKBITS")

    driver = osgeo.gdal.GetDriverByName("GTiff")
    ny, nx = grid.size
    dataset = driver.Create(fnm, nx, ny, 1, gdal_type(grid.values.dtype), co)
    t = grid.transform
    dataset.SetGeoTransform([t[0] - 0.5*(t[2] + t[4]), t[2], -t[4],
                             t[1] + (ny - 0.5)*t[3] + 0.5*t[5], t[5], -t[3]])
    try:
        srs = srs_from_crs(grid.crs)
    except Exception as e:
        sys.stderr.write("Writing GeoTiff failed:\n\t{0}\n".format(e))
        return
    dataset.SetProjection(srs.ExportToWkt())
    band = dataset.GetRasterBand(1)
    band.SetNoDataValue(grid.nodata)
    band.WriteArray(grid.values[::-1])
    band = None
    dataset = None
    return

