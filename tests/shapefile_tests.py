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
from copy import copy

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

        self.points3 = [Point((1, 1, 0)),
                       Point((3, 1, 3)),
                       Point((4, 3, 2)),
                       Point((2, 2, -1))]
        self.line3 = Line([[1,5,2],[5,5,-1],[5,1,3],[3,3,1],[1,1,0]])
        self.polygon3 = Polygon([[1,5,2],[5,5,-1],[5,1,3],[3,3,1],[1,1,0]])
        return

    def test_writepoints(self):
        mp = Multipoint([p.vertex for p in self.points])
        mp.to_shapefile("data/points_shp")
        return

    def test_writeline(self):
        self.line.to_shapefile("data/line_shp")
        return

    def test_writepoly(self):
        self.polygon.to_shapefile("data/polygon_shp")
        return

    def test_writepoints3(self):
        mp = Multipoint([p.vertex for p in self.points3])
        mp.to_shapefile("data/pointsz_shp")
        return

    def test_writeline3(self):
        self.line3.to_shapefile("data/linez_shp")
        return

    def test_writepoly3(self):
        self.polygon3.to_shapefile("data/polygonz_shp")
        return

    def test_write_collection_points(self):
        mp = Multipoint([p.vertex for p in self.points])
        mp0 = copy(mp)
        mp1 = copy(mp.shift((4, 2)))
        mp2 = copy(mp.shift((-2, 3)))
        _shpfuncs.write_shapefile([mp0, mp1, mp2], "data/points_collection")
        return

    def test_write_collection_lines(self):
        line0 = copy(self.line)
        line1 = copy(self.line.shift((4, 2)))
        line2 = copy(self.line.shift((-2, 3)))
        _shpfuncs.write_shapefile([line0, line1, line2], "data/line_collection")
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


