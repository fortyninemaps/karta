"""
Karta is a collection of Python modules for handling vector and raster
geospatial data.
"""

__all__ = ["vector", "raster", "crs", "tile", "errors"]

from .version import __version__

from . import vector
from . import raster
from . import crs
from . import tile
from . import errors

from .vector import (Point, Line, Polygon, Multipoint, Multiline, Multipolygon,
                     from_shape, read_geojson,
                     read_gpx_waypts, read_gpx_tracks,
                     read_shapefile, write_shapefile,
                     geometry)

from .raster import (RegularGrid, SimpleBand, CompressedBand,
                     read_aai, read_geotiff, read_gtiff, from_geotiffs,
                     grid, misc)

# from .tile import tile_bbox, tile_nw_corner, tile_tuple

