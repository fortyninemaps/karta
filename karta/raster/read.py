""" Functions for reading raster data sources as RegularGrid objects """
import numpy as np
from .grid import RegularGrid
from ..crs import Proj4CRS, GeographicalCRS

try:
    from . import _gtiff
    HAS_GDAL = True
except ImportError:
    HAS_GDAL = False

from . import _aai
# from . import _dem

def read_aai(fnm):
    """ Convenience function to open a ESRI ASCII grid and return a RegularGrid
    instance.
    """
    values, aschdr = _aai.aairead(fnm)
    t = {'xllcenter': aschdr['xllcorner'] + 0.5 * aschdr['cellsize'],
         'yllcenter': aschdr['yllcorner'] + 0.5 * aschdr['cellsize'],
         'dx'       : aschdr['cellsize'],
         'dy'       : aschdr['cellsize'],
         'xrot'     : 0.0,
         'yrot'     : 0.0}
    values[values==aschdr['nodata_value']] = np.nan
    return RegularGrid(t, values=values[::-1])

def proj4_isgeodetic(s):
    return ("lonlat" in s) or ("longlat" in s) or \
            ("latlon" in s) or ("latlong" in s)

def read_gtiff(fnm, band=1):
    """ Convenience function to open a GeoTIFF and return a RegularGrid
    instance.

    Parameters
    ----------

    fnm : GeoTiff file path

    band : band to open (default 1)
    """
    if not HAS_GDAL:
        raise NotImplementedError("Right now, loading GeoTiffs requires GDAL.")
    arr, hdr = _gtiff.read(fnm, band)
    t = {'xllcenter'  : hdr['xulcorner'] + 0.5 * (hdr['dx'] + hdr['sx']),
         'yllcenter'  : hdr['yulcorner'] + (hdr['ny'] - 0.5) * hdr['dy'] - 0.5 * hdr['sy'],
         'dx'         : hdr['dx'],
         'dy'         : -hdr['dy'],
         'xrot'       : hdr['sx'],
         'yrot'       : hdr['sy']}

    geodstr = "+a={a} +f={f}".format(a=hdr["srs"]["semimajor"],
                                     f=hdr["srs"]["flattening"])
    if proj4_isgeodetic(hdr["srs"]["proj4"]):
        crs = GeographicalCRS(geodstr, name="Imported GTiff")
    else:
        crs = Proj4CRS(hdr["srs"]["proj4"], geodstr, name="Imported GTiff")
    return RegularGrid(t, values=arr.squeeze()[::-1], crs=crs)

# Aliases for backwards compat.
gtiffread = read_gtiff
aairead = read_aai

