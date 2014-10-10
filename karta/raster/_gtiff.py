""" IO interface to GeoTiffs using GDAL. """

from osgeo import gdal, osr
import osgeo.gdalconst as gc
import numpy as np

gdal.UseExceptions()

def SRS_from_WKT(s):
    """ Return Proj.4 string, semimajor axis, and flattening """
    sr = osr.SpatialReference()
    sr.ImportFromWkt(s)
    return sr

def read(fnm):
    """ Read a GeoTiff file and return a numpy array and a dictionary of header
    information. """
    hdr = dict()
    dataset = gdal.Open(fnm, gc.GA_ReadOnly)

    try:
        hdr["nbands"] = dataset.RasterCount
        hdr["nx"] = dataset.RasterXSize
        hdr["ny"] = dataset.RasterYSize
        arr = np.empty((hdr["nbands"], hdr["ny"], hdr["nx"]))

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

        for i in range(1, hdr["nbands"]+1):
            band = dataset.GetRasterBand(i)
            arr[i-1,:,:] = band.ReadAsArray()

        sr = SRS_from_WKT(dataset.GetProjectionRef())
        hdr["srs"] = {"proj4": sr.ExportToProj4(),
                      "semimajor": sr.GetSemiMajor(),
                      "flattening": sr.GetInvFlattening()}
    finally:
        dataset = None

    return arr, hdr

