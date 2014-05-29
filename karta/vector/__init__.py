"""
Vector data classes and functions.
"""

from . import guppy
from . import xyfile
from . import geojson
from . import gpx
from . import vtk
from . import quadtree

from .guppy import Point, Multipoint, Line, Polygon
from .read import read_geojson, read_geojson_features, read_shapefile
from ._shpfuncs import write_shapefile
from .gpx import GPX
from .xyfile import load_xy, xyz2array_reg, array2xyz

try:
    from . import stats
except ImportError:
    # Probably missing scipy dependency
    pass

__all__ = ["guppy", "xyfile", "geojson", "gpx", "vtk",
           "Point", "Multipoint", "Line", "Polygon",
           "read_geojson", "read_geojson_features", "read_shapefile",
           "write_shapefile",
           "GPX",
           "load_xy", "xyz2array_reg", "array2xyz"]

