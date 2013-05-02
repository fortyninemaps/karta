"""
Classes for handling raster data.
"""

import grid
import aaigrid
import raster
import flow

from grid import RegularGrid, StructuredGrid, aairead
from aaigrid import AAIGrid
from raster import witch_of_agnesi, peaks, pad, slope, aspect, grad, div
from raster import normed_vector_field, fill_sinks

try:
    from cfuncs import streamline2d
except ImportError:
    from streamline import streamline2d

try:
    import fill_sinks    # Cython
except ImportError:
    print("Warning: experimental Cython fill_sinks not loaded")

try:
    import flow          # Has a Scipy dependency
except ImportError:
    print("Warning: raster.flow not loaded - perhaps scipy is missing")

