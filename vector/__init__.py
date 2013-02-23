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
from shp_funcs import read_shapefile, write_shapefile
from gpxparser import GPXParser
from xyfile import load_xy, xyz2array_reg, array2xyz

