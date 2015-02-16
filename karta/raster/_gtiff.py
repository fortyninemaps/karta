""" IO interface to GeoTiffs using GDAL. """

import sys
from osgeo import gdal, osr
import osgeo.gdalconst as gc
import numpy as np

gdal.UseExceptions()

def SRS_from_WKT(s):
    """ Return Proj.4 string, semimajor axis, and flattening """
    sr = osr.SpatialReference()
    sr.ImportFromWkt(s)
    return sr

def numpy_dtype(dt_int):
    """ Return a numpy dtype that matches the band data type. """
    name = gdal.GetDataTypeName(dt_int)
    if name == "UInt16":
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
        return gdal.GDT_UInt16
    elif dtype == np.int16:
        return gdal.GDT_Int16
    elif dtype == np.int32:
        return gdal.GDT_Int32
    elif dtype == np.int32:
        return gdal.GDT_Int32
    elif dtype == np.float32:
        return gdal.GDT_Float32
    elif dtype == np.float64:
        return gdal.GDT_Float64
    elif dtype == np.complex64:
        return gdal.GDT_CFloat64
    else:
        raise TypeError("GDAL equivalent to type {0} unknown".format(dtype))

def read(fnm, band):
    """ Read a GeoTiff file and return a numpy array and a dictionary of header
    information. """
    hdr = dict()
    dataset = gdal.Open(fnm, gc.GA_ReadOnly)

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

        arr = np.empty((hdr["ny"], hdr["nx"]), dtype=numpy_dtype(max_dtype))
        dtype = numpy_dtype(rasterband.DataType)
        arr[:,:] = rasterband.ReadAsArray(buf_obj=np.empty([ny, nx], dtype=dtype))

    finally:
        dataset = None
    return arr, hdr

def srs_from_crs(crs):
    try:
        srs = osr.SpatialReference()
        srs.ImportFromProj4(crs.project.srs)
    except AttributeError:
        raise TypeError("Right now, GDAL projection can only be set from CRS instances backed by proj.4")
    return srs

def write(fnm, grid):
    """ Write a grid-like object to *fnm* """
    driver = gdal.GetDriverByName("GTiff")
    ny, nx = grid.size
    dataset = driver.Create(fnm, nx, ny, 1, gdal_type(grid.values.dtype))
    t = grid.transform
    dataset.SetGeoTransform([t[0] - 0.5*(t[2] + t[4]), t[2], -t[4],
                             t[1] + (ny - 0.5)*t[3] + 0.5*t[5], t[5], -t[3]])
    try:
        srs = srs_from_crs(grid.crs)
    except TypeError as e:
        sys.stderr.write("Writing GeoTiff failed:\n\t{0}\n".format(e))
        return
    dataset.SetProjection(srs.ExportToWkt())
    dataset.GetRasterBand(1).WriteArray(grid.values[::-1])
    dataset = None
    return

