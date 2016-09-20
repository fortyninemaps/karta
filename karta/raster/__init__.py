"""
karta.raster

Classes for grids of raster data
"""

import numpy
import numbers
numbers.Integral.register(numpy.integer)

from . import grid
from . import misc

from .grid import RegularGrid, merge, gridpoints, mask_poly
from .band import SimpleBand, CompressedBand
from .read import read_aai, read_geotiff, read_gtiff, from_geotiffs
from .misc import (normed_potential_vectors,
                   slope, aspect, gradient, divergence, hillshade)

__all__ = ["grid", "misc", "RegularGrid",
           "read_aai", "read_geotiff", "read_gtiff", "from_geotiffs",
           "slope", "aspect", "gradient", "divergence", "hillshade",
           "normed_potential_vectors"]

