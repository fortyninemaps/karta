"""
Karta is a collection of Python modules for handling vector and raster
geospatial data.
"""

import pkgutil
__path__ = pkgutil.extend_path(__path__, __name__)

__all__ = ["vector", "raster", "crs", "errors"]
from .version import __version__

from . import vector
from . import raster
from . import crs
from . import errors
from .vector import *
from .raster import *

