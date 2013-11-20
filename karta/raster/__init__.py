"""
Classes for handling raster data.
"""

import grid
import aaigrid
import raster
import flow

from grid import RegularGrid, WarpedGrid, aairead
from aaigrid import AAIGrid
from raster import witch_of_agnesi, peaks, pad, slope, aspect, grad, div
from raster import normed_vector_field

try:
    from crfuncs import streamline2d
except ImportError:
    from streamline import streamline2d

try:
    import cfill_sinks as fill_sinks
except ImportError:
    import fill_sinks

try:
    import flow          # Has a Scipy dependency
except ImportError:
    pass

__all__ = ["grid", "aaigrid", "raster", "flow",
           "RegularGrid", "WarpedGrid", "aairead",
           "AAIGrid",
           "pad", "slope", "aspect", "grad", "div", "normed_vector_field",
           "streamline2d",
           "fill_sinks"]

