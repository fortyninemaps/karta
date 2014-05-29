""" IO interface to GeoTiffs using GDAL. """

import gdal
import gdalconst as gc
import numpy as np

gdal.UseExceptions()

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
            hdr["dy"] = -transform[5]
            hdr["xllcorner"] = transform[0]
            hdr["yllcorner"] = transform[3] - hdr["ny"] * hdr["dy"]
        else:
            raise AttributeError("No GeoTransform in GDAL Dataset")

        for i in range(1, hdr["nbands"]+1):
            band = dataset.GetRasterBand(i)
            arr[i-1,:,:] = band.ReadAsArray()

        projstr = dataset.GetProjection()
        hdr["crs"] = projstr

    finally:
        dataset = None

    return arr, hdr

