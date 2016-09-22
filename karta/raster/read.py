""" Functions for reading raster data sources as RegularGrid objects """
import numpy as np
from .grid import RegularGrid
from ..crs import ProjectedCRS, GeographicalCRS
from . import _gdal
from . import _aai
from ..errors import GridError

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

def read_gtiff(fnm, in_memory=True, ibands=_gdal.ALL, **kw):
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
        will have bands of type karta.raster._gdal.GdalFileBand
    """
    bands, hdr = _gdal.read(fnm, in_memory, ibands, **kw)

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

def from_geotiffs(*fnms, **kw):
    """ Read multiple GeoTiff files as bands within a single grid. Reads the
    first band within each GeoTiff, and checks that grids share the same
    spatial transform and coordinate system.

    Parameters
    ----------
    fnm : list of str
        GeoTiff file paths
    in_memory : bool
        if True (default), read entire dataset into memory
    bandclass : Band class, optional
        class of band used by returned grid (default karta.band.CompressedBand)
        if in_memory is True, this parameter is ignored and the returned grid
        will have bands of type karta.raster._gdal.GdalFileBand
    """
    in_memory = kw.get("in_memory", True)
    if len(fnms) == 0:
        raise ValueError("at least one filename must be provided")
    bands = []
    hdrs = []
    for i, fnm in enumerate(fnms):
        _b, _h = _gdal.read(fnm, in_memory, 1, **kw)
        if (i != 0) and (len(hdrs) != 0):
            for k, v in _h.items():
                if (isinstance(v, float)
                        and np.isnan(v) and not np.isnan(hdrs[-1].get(k, 0))):
                    if v != hdrs[-1].get(k, np.nan):
                        raise GridError("geotransform in {0} not equivalent to "
                                        "a previous GeoTiff".format(fnm))
        bands.append(_b[0])
        hdrs.append(_h)

    hdr = hdrs[0]
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
read_geotiff = read_gtiff
aairead = read_aai

