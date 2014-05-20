import unittest
import os
import math
import numpy as np
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
from test_helper import md5sum, md5sum_file, TESTDATA, TESTDIR

import datetime
import numbers

from karta.vector import _shpfuncs
import karta.vector as vector
import karta.crs as crs
from karta.vector.guppy import Point, Multipoint, Line, Polygon

class TestShapefile(unittest.TestCase):

    def setUp(self):
        self.points = [Point((1, 1)),
                       Point((3, 1)),
                       Point((4, 3)),
                       Point((2, 2))]
        self.line = Line([[1,5],[5,5],[5,1],[3,3],[1,1]])
        self.polygon = Polygon([[1,5],[5,5],[5,1],[3,3],[1,1]])
        return

    def test_writeline(self):
        self.line.to_shapefile("data/line_shp")
        return

    def test_dbase_type(self):
        self.assertEqual(_shpfuncs.property_field_type(1.0), "O")
        self.assertEqual(_shpfuncs.property_field_type(1), "I")
        #self.assertEqual(_shpfuncs.property_field_type(True), "L")
        #self.assertEqual(_shpfuncs.property_field_type(False), "L")
        self.assertEqual(_shpfuncs.property_field_type("pale ale"), "C")
        self.assertEqual(_shpfuncs.property_field_type(datetime.date(1986, 8, 17)), "D")
        self.assertEqual(_shpfuncs.property_field_type(datetime.datetime(2013, 5, 4, 20, 40, 21)), "@")
        return

if __name__ == "__main__":
    unittest.main()


