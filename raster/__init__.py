"""
Classes for handling raster data.
"""

import aaigrid
from raster import *
import flow

try:
    from cfuncs import streamline2d
except ImportError:
    from streamline import streamline2d
