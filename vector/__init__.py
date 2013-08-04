"""
Vector data classes and functions.
"""

import guppy
import stats
import xyfile
import geojson
import gpx
import vtk
import shapefile
import shp_funcs

from guppy import Point, Multipoint, Line, Polygon
from read import read_geojson, read_geojson_features, read_shapefile
from shp_funcs import write_shapefile
from gpx import GPX
from xyfile import load_xy, xyz2array_reg, array2xyz

try:
    import stats
except ImportError:
    # Probably missing scipy dependency
    pass

__all__ = ["guppy", "xyfile", "geojson", "gpx", "vtk", "shapefile", "shp_funcs",
           "Point", "Multipoint", "Line", "Polygon",
           "read_geojson", "read_geojson_features", "read_shapefile",
           "write_shapefile",
           "GPX",
           "load_xy", "xyz2array_reg", "array2xyz"]

