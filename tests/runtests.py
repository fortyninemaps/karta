#! /usr/bin/env python
""" Run all tests """

import os
import shutil
import unittest

from test_helper import TMPDATA

# Core geometry tests
from crs_tests import *
from coordstring_tests import *
from new_geometry_tests import *
from geometry_tests import *
from dateline_tests import *
from table_tests import *

from grid_tests import *
from band_tests import *
from raster_misc_tests import *

# Vector operations
from vector_predicate_tests import *
from quadtree_tests import *
from rtree_tests import *
#from redblack_tests import *

# Vector IO
from shapefile_tests import *
from geojson_tests import *
from gpx_tests import *
from geointerface_tests import *
from tile_tests import *

# Raster IO
from crfuncs_tests import *
from geotiff_tests import *
from dem_tests import *
#from aai_tests import *

if __name__ == "__main__":

    if not os.path.isdir(TMPDATA):
        os.mkdir(TMPDATA)

    try:
        unittest.main()
    finally:
        if os.path.isdir(TMPDATA):
            shutil.rmtree(TMPDATA)
