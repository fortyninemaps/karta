"""
Classes for handling raster data.
"""

from . import grid
from . import aaigrid
from . import raster

from .grid import RegularGrid, WarpedGrid, aairead, gtiffread
from .aaigrid import AAIGrid
from .raster import witch_of_agnesi, peaks, pad, slope, aspect, grad, div
from .raster import normed_vector_field

try:
    from .crfuncs import streamline2d
except ImportError:
    from .streamline import streamline2d

#try:
#    from . import cfill_sinks as fill_sinks
#except ImportError:
#    from . import fill_sinks

#try:
#    from . import flow          # Has a Scipy dependency
#except ImportError:
#    pass

__all__ = ["grid", "aaigrid", "raster",
           "RegularGrid", "WarpedGrid", "aairead", "gtiffread",
           "AAIGrid",
           "pad", "slope", "aspect", "grad", "div", "normed_vector_field",
           "streamline2d",
           "fill_sinks"]

