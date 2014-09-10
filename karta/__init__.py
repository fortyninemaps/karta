"""
Karta is a collection of Python modules for handling vector and raster
geospatial data.
"""

__version__ = 0.4
__all__ = ["vector", "raster", "crs"]

from . import vector
from . import raster
from . import crs as _crs
from .crs import crsreg as crs
from .vector import *
from .raster import *
from .crs import *

