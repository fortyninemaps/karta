#! /usr/bin/env python
""" Run all tests """

import unittest
import os
from test_helper import TESTDIR
from vector_tests import *
from vector_io_tests import *
from raster_tests import *
from shapefile_tests import *

if __name__ == "__main__":

    TMPDATA = os.path.join(TESTDIR, "data")
    if not os.path.isdir(TMPDATA):
        print("creating {0}".format(TMPDATA))
        os.mkdir(TMPDATA)

    unittest.main()

