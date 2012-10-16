"""
Classes for handling raster data.
"""

import aaigrid_driver
from raster import *
import flow

try:
    from cfuncs import streamline2d
except ImportError:
    from streamline import streamline2d
