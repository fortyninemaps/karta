"""
Karta is a collection of Python modules for handling vector and raster
geospatial data.
"""

import pkgutil
__path__ = pkgutil.extend_path(__path__, __name__)

__all__ = ["vector", "raster", "crs", "tile", "errors"]
from .version import __version__

from . import vector
from . import raster
from . import crs
from . import tile
from . import errors
from .vector import *
from .raster import *
# from .tile import tile_bbox, tile_nw_corner, tile_tuple

