""" Functions for reading raster data sources as RegularGrid objects """
import numpy as np
from .grid import RegularGrid
from ..crs import Proj4CRS, GeographicalCRS
from . import _gtiff
from . import _aai
# from . import _dem

def read_aai(fnm):
    """ Convenience function to open a ESRI ASCII grid and return a RegularGrid
    instance.
    """
    values, aschdr = _aai.aairead(fnm)
    t = {'xllcorner': aschdr['xllcorner'],
         'yllcorner': aschdr['yllcorner'],
         'dx'       : aschdr['cellsize'],
         'dy'       : aschdr['cellsize'],
         'xrot'     : 0.0,
         'yrot'     : 0.0}
    values[values==aschdr['nodata_value']] = np.nan
    return RegularGrid(t, values=values[::-1])

def proj4_isgeodetic(s):
    return ("lonlat" in s) or ("longlat" in s) or \
            ("latlon" in s) or ("latlong" in s)

def read_gtiff(fnm, band=1, in_memory=True):
    """ Convenience function to open a GeoTIFF and return a RegularGrid
    instance.

    Parameters
    ----------

    fnm : GeoTiff file path

    band : band to open (default 1)

    in_memory : if True (default), read entire dataset into memory
    """
    arr, hdr = _gtiff.read(fnm, band, in_memory)
    if isinstance(arr, np.ndarray):
        arr = arr.squeeze()[::-1]

    t = {'xllcorner'  : hdr['xulcorner'] - hdr['ny'] * hdr['sx'],
         'yllcorner'  : hdr['yulcorner'] + hdr['ny'] * hdr['dy'],
         'dx'         : hdr['dx'],
         'dy'         : -hdr['dy'],
         'xrot'       : hdr['sx'],
         'yrot'       : -hdr['sy']}

    geodstr = "+a={a} +f={f}".format(a=hdr["srs"]["semimajor"],
                                     f=hdr["srs"]["flattening"])
    if proj4_isgeodetic(hdr["srs"]["proj4"]):
        crs = GeographicalCRS(geodstr, name="Imported GTiff")
    else:
        crs = Proj4CRS(hdr["srs"]["proj4"], geodstr, name="Imported GTiff")
    return RegularGrid(t, values=arr, crs=crs, nodata_value=hdr["nodata"])

# Aliases for backwards compat.
gtiffread = read_gtiff
aairead = read_aai

