"""
Karta is a collection of python modules for handling vector and raster
geospatial data.
"""

import vector
import vector.stats
import vector.guppy as guppy
import vector.shapefile as shapefile
import vector.xyfile as xyfile
import vector.geojson
import vector.vtk

import raster
import raster.raster
import raster.streamline
import raster.aaigrid as aaigrid

try:
    import raster.fill_sinks    # Cython
except ImportError:
    print("Warning: experimental Cython fill_sinks not loaded")
try:
    import raster.flow          # Has a Scipy dependency
except ImportError:
    print("Warning: raster.flow not loaded - perhaps scipy is missing")

__all__ = ["vector", "raster"]

