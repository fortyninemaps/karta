"""
karta.raster

Classes for grids of raster data
"""

import numpy
import numbers
numbers.Integral.register(numpy.integer)

from . import grid
from . import misc

from .grid import RegularGrid, WarpedGrid, merge, gridpoints, mask_poly
from .band import SimpleBand, CompressedBand
from .read import read_aai, read_gtiff, from_geotiffs, aairead, gtiffread
from .misc import (witch_of_agnesi, pad, normed_potential_vectors,
                   slope, aspect, gradient, divergence, hillshade)

__all__ = ["grid", "misc",
           "RegularGrid", "WarpedGrid",
           "aairead", "gtiffread", "read_aai", "read_gtiff", "from_geotiffs",
           "slope", "aspect", "gradient", "divergence", "hillshade",
           "normed_potential_vectors"]

