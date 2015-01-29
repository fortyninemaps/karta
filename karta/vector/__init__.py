"""
Vector data classes and functions.
"""

from . import geometry
from . import xyfile
from . import geojson
from . import gpx
from . import vtk
from . import quadtree

from .geometry import Point, Multipoint, Line, Polygon
from .read import read_geojson, read_shapefile
from .shp import write_shapefile
from .gpx import GPX
from .xyfile import load_xy, xyz2array_reg, array2xyz

__all__ = ["geometry", "xyfile", "geojson", "gpx", "vtk",
           "Point", "Multipoint", "Line", "Polygon",
           "read_geojson", "read_shapefile",
           "write_shapefile",
           "GPX",
           "load_xy", "xyz2array_reg", "array2xyz"]

