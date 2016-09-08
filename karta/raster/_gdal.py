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

    def __getitem__(self, idx):
        ny, nx = self.size

        if isinstance(idx, tuple):
            iidx, jidx = idx
        else:
            iidx = idx
            jidx = slice(0, nx, 1)

        # Calculate start, end, step ranges for rows and columns
        if isinstance(iidx, int):
            ystart = iidx
            yend = iidx+1
            ystep = 1
        else:
            ystart, yend, ystep = iidx.indices(ny)

        if isinstance(jidx, int):
            xstart = jidx
            xend = jidx+1
            xstep = 1
        else:
            xstart, xend, xstep = jidx.indices(nx)

        # Extract ...
        if abs(yend-ystart) == 1:
            # ... a row vector
            x0 = min(xstart, xend)
            y0 = min(ny-ystart, ny-yend)
            ret = self.gdalband.ReadAsArray(x0, y0, abs(xend-xstart), 1)[0,::xstep]

            if abs(xend-xstart) == 1:
                ret = ret[0]

        elif abs(xend-xstart) == 1:
            # ... a column vector
            x0 = min(xstart, xend)
            y0 = min(ny-ystart, ny-yend)
            ret = self.gdalband.ReadAsArray(x0, y0, 1, abs(yend-ystart))[::ystep].ravel()

        elif (abs(xstep) == 1) and (abs(ystep) == 1):
            # ... contiguous blocks (fast path)
            x0 = min(xstart, xend)
            y0 = min(ny-ystart, ny-yend)
            ret = self.gdalband.ReadAsArray(x0, y0, abs(xend-xstart), abs(yend-ystart))
            if xstep < 0:
                ret = ret[:,::-1]

            if ystep > 0:       # Flip from GDAL convention to karta convention
                ret = ret[::-1]

        else:
            # Sparse extraction
            t = self.gdalband.DataType
            values = (self.gdalband.ReadRaster(x, ny-y-1, 1, 1, 1, 1, t)
                        for y in range(ystart, yend, ystep)
                        for x in range(xstart, xend, xstep))

            # compute size of output
            outnx = ceil(float(abs(xend - xstart)) / abs(xstep))
            outny = ceil(float(abs(yend - ystart)) / abs(ystep))

            pyt = pytype(t)
            values_py = [struct.unpack(pyt, a) for a in values]
            ret = np.array(values_py).reshape([outny, outnx])
        return ret

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
    sr.ImportFromWkt(s)
    return sr

def pytype(dt_int):
    """ Return a fmt character based on a gdal type integer. """
    if dt_int == 1:
        return "B"
    elif dt_int == 3:
        return "h"
    elif dt_int == 2:
        return "H"
    elif dt_int == 5:
        return "i"
    elif dt_int == 4:
        return "I"
    #elif dt_int == 5:
    #    return "l"
    #elif dt_int == 4:
    #    return "L"
    elif dt_int == 6:
        return "f"
    elif dt_int == 7:
        return "d"
    else:
        raise KeyError("No Python type to coerce from GDT_Int %d" % dt_int)

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
    if dtype == np.uint8:
        return osgeo.gdal.GDT_Byte
    elif dtype == np.uint16:
        return osgeo.gdal.GDT_UInt16
    elif dtype == np.int8:
        return osgeo.gdal.GDT_Byte      # transform -127 -- 127 to 0 -- 255
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
        hdr["srs"] = {"proj4": sr.ExportToProj4(),
                      "semimajor": sr.GetSemiMajor(),
                      "flattening": sr.GetInvFlattening(),
                      "name": sr.GetAttrValue('PROJCS')}

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
                bands[i][:,:] = _arr.squeeze()[::-1]
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
    try:
        srs = srs_from_crs(grid.crs)
    except Exception as e:
        sys.stderr.write("Writing GeoTiff failed:\n\t{0}\n".format(e))
        return
    dataset.SetProjection(srs.ExportToWkt())
    for i, _ in enumerate(grid.bands):
        band = dataset.GetRasterBand(i+1)
        band.SetNoDataValue(grid.nodata)
        band.WriteArray(grid.bands[i][::-1])
    band = None
    dataset = None
    return

