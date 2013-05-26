"""
Vector data classes and functions.
"""

import guppy
import stats
import xyfile
import geojson
import vtk
import shapefile
import shp_funcs

from guppy import Point, Multipoint, Line, Polygon
from read import read_geojson, read_geojson_features, read_shapefile
from shp_funcs import write_shapefile
from gpxparser import GPXParser
from xyfile import load_xy, xyz2array_reg, array2xyz

try:
    import stats
except ImportError:
    # Probably missing scipy dependency
    pass

