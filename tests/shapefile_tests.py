import unittest
import os
import math
import numpy as np
from test_helper import md5sum, md5sum_file, TESTDATA, TESTDIR

import datetime
import numbers
from copy import copy

from karta.vector import _shpfuncs, read_shapefile
from karta.vector.guppy import Point, Multipoint, Line, Polygon

class TestShapefile(unittest.TestCase):

    def setUp(self):
        self.points = [Point((1, 1), data={"species": "T. officianale"}),
                       Point((3, 1), data={"species": "C. tectorum"}),
                       Point((4, 3), data={"species": "M. alba"}),
                       Point((2, 2), data={"species": "V. cracca"})]
        self.multipoint = Multipoint([(1,1), (3,1), (4,3), (2,2)],
                                     data={"species": ["T. officianale", "C. tectorum",
                                                       "M. alba", "V. cracca"]})
        self.line = Line([(1.0,5.0),(5.0,5.0),(5.0,1.0),(3.0,3.0),(1.0,1.0)])
        self.polygon = Polygon([(1.0,5.0),(5.0,5.0),(5.0,1.0),(3.0,3.0),(1.0,1.0)])

        self.points3 = [Point((1, 1, 0)),
                       Point((3, 1, 3)),
                       Point((4, 3, 2)),
                       Point((2, 2, -1))]
        self.line3 = Line([(1,5,2),(5,5,-1),(5,1,3),(3,3,1),(1,1,0)])
        self.polygon3 = Polygon([(1,5,2),(5,5,-1),(5,1,3),(3,3,1),(1,1,0)])

        exists = os.path.exists
        join = os.path.join
        testfiles = ["points.shp", "line.shp", "polygon.shp"]
        if any(not exists(join(TESTDATA, "shapefiles/", fnm)) for fnm in testfiles):
            self.saveTestData()
        return

    def saveTestData(self):
        testfiles = [(self.multipoint, "points"),
                     (self.line, "line"),
                     (self.polygon, "polygon")]
        for (geom, fnm) in testfiles:
            geom.to_shapefile(os.path.join(TESTDATA, "shapefiles", fnm))
        return

    def assertGeomEqual(self, this, that):
        self.assertTrue(np.all(this.get_vertices() == that.get_vertices()))
        self.assertEqual(this._crs, that._crs)
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
        self.assertEqual(_shpfuncs.property_field_type(np.float32(1.0)), "O")
        self.assertEqual(_shpfuncs.property_field_type(np.int16(1)), "I")
        #self.assertEqual(_shpfuncs.property_field_type(True), "L")
        #self.assertEqual(_shpfuncs.property_field_type(False), "L")
        self.assertEqual(_shpfuncs.property_field_type("pale ale"), "C")
        self.assertEqual(_shpfuncs.property_field_type(datetime.date(1986, 8, 17)), "D")
        self.assertEqual(_shpfuncs.property_field_type(datetime.datetime(2013, 5, 4, 20, 40, 21)), "@")
        return

    def test_read_points(self):
        shps = read_shapefile(os.path.join(TESTDATA, "newp"))
        mp = shps[0]
        self.assertEqual(mp.vertices,
                         [(-14.612, 80.50906666666667), (-14.612, 80.50906666666667),
                          (-14.612, 80.50906666666667), (-13.744733333333333, 80.28181666666667),
                          (-13.744733333333333, 80.28181666666667), (-13.744733333333333, 80.28181666666667),
                          (-11.002583333333334, 80.32173333333333), (-11.002583333333334, 80.32173333333333),
                          (-11.002583333333334, 80.32173333333333), (-11.07225, 80.56316666666666),
                          (-11.07225, 80.56316666666666)])
        self.assertEqual(mp.data["meterno"], ['IMS1/1', 'IMS2/1', '5952/2',
                                              'IMS4/1', '5953/2', '1963/13',
                                              'IMS5/1', '5213/A', '2121/13',
                                              'IMS3/1', '3613/2'])
        self.assertEqual(mp.data["depth_m"], [73, 143, 247, 86, 147, 250, 74,
                                              142, 235, 150, 248])
        return

    def test_read_multipoint_attributes(self):
        mp = read_shapefile(os.path.join(TESTDATA, "shapefiles", "points"))
        self.assertEqual(mp[0].data["species"], self.multipoint.data["species"])
        return

    def test_read_line(self):
        line = read_shapefile(os.path.join(TESTDATA, "shapefiles", "line"))[0]
        self.assertGeomEqual(line, self.line)
        return

    def test_read_polygon(self):
        polygon = read_shapefile(os.path.join(TESTDATA, "shapefiles", "polygon"))[0]
        self.assertGeomEqual(polygon, self.polygon)
        return

if __name__ == "__main__":
    unittest.main()


