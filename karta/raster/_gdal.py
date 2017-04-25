""" IO interface to GeoTiffs using GDAL. """

import struct
import sys
from math import ceil
import numpy as np
from .band import CompressedBand
from .. import errors

import osgeo.gdal
import osgeo.osr
import osgeo.gdalconst as gc
osgeo.gdal.UseExceptions()

ALL = -1

class GdalFileBand(object):
    """ Read-only raster Band interface the reads data from a disk-bound
    datasource.
    """
    def __init__(self, band, dataset):
        """
        Parameters
        ----------
        band : osgeo.gdal.Band
        dataset : osgeo.gdal.Dataset
        """
        self.gdalband = band
        self.dataset = dataset
        return

    def __del__(self):
        self.dataset = None
        self.gdalband = None

    def getblock(self, yoff, xoff, ny, nx):
        # Note that GDAL uses the alternative x,y convention
        grid_ny, grid_nx = self.size
        chunk = self.gdalband.ReadAsArray(xoff, grid_ny - yoff - ny, nx, ny)
        if chunk is None:
            raise IOError("failure reading slice from GDAL backend")
        return chunk[::-1]

    def setblock(self, yoff, xoff, array):
        raise NotImplementedError()

    def __iter__(self):
        for i in range(self.dataset.RasterYSize):
            yield self[i,:]

    @property
    def size(self):
        return (self.dataset.RasterYSize, self.dataset.RasterXSize)

    @property
    def dtype(self):
        return np.dtype(numpy_dtype(self.gdalband.DataType))

def SRS_from_WKT(s):
    """ Return Proj.4 string, semimajor axis, and flattening """
    sr = osgeo.osr.SpatialReference()
    try:
        sr.ImportFromWkt(s)
    except RuntimeError:
        sr = None
    return sr

NUMPY_DTYPE = {"Byte": np.uint8,
               "UInt16": np.uint16,
               "Int16": np.int16,
               "UInt32": np.uint32,
               "Int32": np.int32,
               "Float32": np.float32,
               "Float64": np.float64,
               "CInt16": np.complex64,
               "CInt32": np.complex64,
               "CFloat32": np.complex64,
               "CFloat64": np.complex64}

def numpy_dtype(dt_int):
    """ Return a numpy dtype that matches the band data type. """
    name = osgeo.gdal.GetDataTypeName(dt_int)
    if name in NUMPY_DTYPE:
        return NUMPY_DTYPE[name]
    else:
        raise TypeError("GDAL data type {0} unknown to karta".format(dt_int))

def gdal_type(dtype):
    """ Return a GDAL type that most closely matches numpy dtype

    Notes
    -----
    Returns GDT_Int32 for np.int64, which may result in overflow.
    """
    if dtype == np.uint8:
        return osgeo.gdal.GDT_Byte
    elif dtype == np.uint16:
        return osgeo.gdal.GDT_UInt16
    elif dtype == np.int8:
        return osgeo.gdal.GDT_Byte      # transform -127 -- 127 to 0 -- 255
    elif dtype == np.int16:
        return osgeo.gdal.GDT_Int16
    elif (dtype == np.int32) or (dtype == np.int64):
        return osgeo.gdal.GDT_Int32
    elif dtype == np.float32:
        return osgeo.gdal.GDT_Float32
    elif dtype == np.float64:
        return osgeo.gdal.GDT_Float64
    elif dtype == np.complex64:
        return osgeo.gdal.GDT_CFloat64
    else:
        raise TypeError("GDAL equivalent to type {0} unknown".format(dtype))

def read(fnm, in_memory, ibands=ALL, bandclass=CompressedBand):
    """ Read a GeoTiff file and return a numpy array and a dictionary of header
    information.

    Parameters
    ----------
    fnm : str
        input datasource
    in_memory : boolean
        indicates whether array should be read fully into memory
    ibands : int or list of ints
        band number (1...)
    bandclass : karta.raster.band class
        if *in_memory* is `False`, use this class for band storage

    Returns an band object and a dictionary of metadata
    """
    hdr = dict()
    dataset = osgeo.gdal.Open(fnm, gc.GA_ReadOnly)

    if ibands == ALL:
        ibands = list(range(1, dataset.RasterCount+1))
    elif not hasattr(ibands, "__iter__"):
        ibands = [ibands]

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
        if sr is not None:
            hdr["srs"] = {"proj4": sr.ExportToProj4(),
                          "semimajor": sr.GetSemiMajor(),
                          "flattening": sr.GetInvFlattening(),
                          "name": sr.GetAttrValue('PROJCS')}
        else:
            hdr["srs"] = {"proj4": "",
                          "semimajor": 6370997.0,
                          "flattening": 1.0 / 298.257223563,
                          "name": "NA"}

        max_dtype = 0
        rasterbands = [dataset.GetRasterBand(i) for i in ibands]
        hdr["nodata"] = rasterbands[0].GetNoDataValue()
        nx = rasterbands[0].XSize
        ny = rasterbands[0].YSize
        if rasterbands[0].DataType > max_dtype:
            max_dtype = rasterbands[0].DataType

        if in_memory:
            dtype = numpy_dtype(rasterbands[0].DataType)
            bands = [bandclass((ny, nx), dtype) for _ in ibands]
            for i, rb in enumerate(rasterbands):
                _arr = rb.ReadAsArray(buf_obj=np.empty([ny, nx], dtype=dtype))
                if _arr is None:
                    raise IOError("error reading GDAL band {}".format(i+1))
                bands[i].setblock(0, 0, _arr.squeeze()[::-1])
        else:
            bands = [GdalFileBand(rb, dataset) for rb in rasterbands]

    finally:
        if in_memory:
            dataset = None
    return bands, hdr

def srs_from_crs(crs):
    srs = osgeo.osr.SpatialReference()
    # SpatialReference can't parse 'lonlat'
    proj4 = crs.get_proj4().replace("lonlat", "latlong")
    srs.ImportFromProj4(proj4)
    return srs

def write(fnm, grid, compress=None, tiled=False, **kw):
    """ Write a grid-like object to *fnm* """
    co = []
    if compress == "LZW":
        co.append("COMPRESS=LZW")
    elif compress == "PACKBITS":
        co.append("COMPRESS=PACKBITS")
    elif compress == "DEFLATE":
        co.append("COMPRESS=DEFLATE")
    elif compress == "LZMA":
        co.append("COMPRESS=LZMA")

    if tiled:
        co.append("TILED=YES")

    for k, v in kw.items():
        co.append("{0}={1}".format(k,v))

    driver = osgeo.gdal.GetDriverByName("GTiff")
    ny, nx = grid.size
    dataset = driver.Create(fnm, nx, ny, len(grid.bands), gdal_type(grid.values.dtype), co)
    t = grid.transform
    dataset.SetGeoTransform([t[0] + ny*t[4], t[2], -t[4],
                             t[1] + ny*t[3], t[5], -t[3]])
    srs = srs_from_crs(grid.crs)
    dataset.SetProjection(srs.ExportToWkt())
    for i, _ in enumerate(grid.bands):
        band = dataset.GetRasterBand(i+1)
        band.SetNoDataValue(grid.nodata)

        tmp = grid[:,:,i]
        band.WriteArray(tmp)
        grid[:,:,i] = tmp[::-1]
    band = None
    dataset = None
    return

