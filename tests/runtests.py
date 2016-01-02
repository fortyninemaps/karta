#! /usr/bin/env python
""" Run all tests """

import unittest
import os
from test_helper import TESTDIR

# Core geometry tests
from crs_tests import *
from geometry_init_tests import *
from geometry_tests import *
from quadtree_tests import *
from rtree_tests import *
from metadata_tests import *
from raster_tests import *

# Vector IO
from shapefile_tests import *
from geojson_tests import *
from misc_io_tests import *

# Raster IO
from geotiff_tests import *

if __name__ == "__main__":

    TMPDATA = os.path.join(TESTDIR, "data")
    if not os.path.isdir(TMPDATA):
        print("creating {0}".format(TMPDATA))
        os.mkdir(TMPDATA)

    unittest.main()
