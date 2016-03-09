"""
Classes for handling raster data.
"""

from . import grid
from . import misc

from .grid import RegularGrid, WarpedGrid, merge, gridpoints, mask_poly
from .read import read_aai, read_gtiff, aairead, gtiffread
from .misc import (witch_of_agnesi, peaks, pad, normed_potential_vectors,
                   slope, aspect, gradient, divergence, hillshade)
from .crfuncs import streamline2d

__all__ = ["grid", "misc",
           "RegularGrid", "WarpedGrid",
           "aairead", "gtiffread", "read_aai", "read_gtiff",
           "slope", "aspect", "gradient", "divergence", "hillshade",
           "normed_potential_vectors", "streamline2d"]

