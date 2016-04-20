""" Functions for reading raster data sources as RegularGrid objects """
import numpy as np
from .grid import RegularGrid
from ..crs import ProjectedCRS, GeographicalCRS
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
         'sx'       : 0.0,
         'sy'       : 0.0}
    values[values==aschdr['nodata_value']] = np.nan
    return RegularGrid(t, values=values[::-1])

def proj4_isgeodetic(s):
    return ("lonlat" in s) or ("longlat" in s) or \
            ("latlon" in s) or ("latlong" in s)

def read_gtiff(fnm, in_memory=True, ibands=_gtiff.ALL, **kw):
    """ Convenience function to open a GeoTIFF and return a RegularGrid
    instance.

    Parameters
    ----------
    fnm : str
        GeoTiff file path
    in_memory : bool
        if True (default), read entire dataset into memory
    ibands : int | list of ints, optional
        band(s) to open (default all)
    bandclass : Band class, optional
        class of band used by returned grid (default karta.band.CompressedBand)
        if in_memory is True, this parameter is ignored and the returned grid
        will have bands of type karta.raster._gtiff.GdalFileBand
    """
    bands, hdr = _gtiff.read(fnm, in_memory, ibands, **kw)

    t = {'xllcorner': hdr['xulcorner'] - hdr['ny'] * hdr['sx'],
         'yllcorner': hdr['yulcorner'] + hdr['ny'] * hdr['dy'],
         'dx'       : hdr['dx'],
         'dy'       : -hdr['dy'],
         'sx'       : hdr['sx'],
         'sy'       : -hdr['sy']}

    if proj4_isgeodetic(hdr["srs"]["proj4"]):
        geodstr = "+a={a} +f={f}".format(a=hdr["srs"]["semimajor"],
                                         f=hdr["srs"]["flattening"])
        crs = GeographicalCRS(geodstr, name=hdr["srs"]["name"])
    else:
        crs = ProjectedCRS(hdr["srs"]["proj4"], name=hdr["srs"]["name"])
    return RegularGrid(t, bands=bands, crs=crs, nodata_value=hdr["nodata"])

# Aliases for backwards compat.
gtiffread = read_gtiff
aairead = read_aai

