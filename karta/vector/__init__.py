"""
karta.vector

Classes and functions for vector data
"""

import numpy
import numbers
numbers.Integral.register(numpy.integer)

from . import geometry
from . import quadtree
from . import rtree
from . import table
from . import _geojson as geojson
from . import _gpx as gpx
from . import _shp as shp
#from . import redblack

from .geometry import (Geometry, Point, Line, Polygon,
                       Multipoint, Multiline, Multipolygon,
                       multipart_from_singleparts, merge_multiparts)
from .read import (from_shape, read_geojson, read_shapefile,
                   read_gpx_waypts, read_gpx_tracks)
from .table import Table
from ._shp import write_shapefile

__all__ = ["geometry", "geojson", "gpx",
           "Point", "Line", "Polygon",
           "Multipoint", "Multiline", "Multipolygon",
           "multipart_from_singleparts",
           "Table",
           "read_geojson", "read_shapefile", "write_shapefile"]

