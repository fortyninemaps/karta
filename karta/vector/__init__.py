"""
karta.vector

Classes and functions for vector data
"""

import numpy
import numbers
numbers.Integral.register(numpy.integer)

from . import geometry
from . import table
from . import xyfile
from . import geojson
from . import gpx
from . import quadtree
from . import rtree
#from . import redblack

from .geometry import (Geometry, Point, Line, Polygon,
                       Multipoint, Multiline, Multipolygon,
                       multipart_from_singleparts)
from .read import (from_shape, read_geojson, read_shapefile,
                   read_gpx_waypts, read_gpx_tracks)
from .shp import write_shapefile
from .table import Table

__all__ = ["geometry", "xyfile", "geojson", "gpx",
           "Point", "Line", "Polygon",
           "Multipoint", "Multiline", "Multipolygon",
           "multipart_from_singleparts",
           "Table",
           "read_geojson", "read_shapefile", "write_shapefile"]

