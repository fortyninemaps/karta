import unittest
import os
from os.path import exists, join
import numpy as np
from test_helper import TESTDIR, TESTDATA, TMPDATA
import datetime
from copy import copy

from karta.vector import shp, read_shapefile
from karta.vector.geometry import Point, Multipoint, Line, Polygon
from karta.crs import LonLatWGS84

class TestShapefile(unittest.TestCase):

    def setUp(self):
        self.points = [Point((1, 1), data={"species": "T. officianale"}, crs=LonLatWGS84),
                       Point((3, 1), data={"species": "C. tectorum"}, crs=LonLatWGS84),
                       Point((4, 3), data={"species": "M. alba"}, crs=LonLatWGS84),
                       Point((2, 2), data={"species": "V. cracca"}, crs=LonLatWGS84)]

        self.multipoint = Multipoint([(1,1), (3,1), (4,3), (2,2)],
                                     data={"species": ["T. officianale", "C. tectorum",
                                                       "M. alba", "V. cracca"]},
                                     crs=LonLatWGS84)

        self.line = Line([(1.0,5.0),(5.0,5.0),(5.0,1.0),(3.0,3.0),(1.0,1.0)],
                         properties={"geom_id": 27, "name": "test line"},
                         crs=LonLatWGS84)

        self.polygon = Polygon([(1.0,5.0),(5.0,5.0),(5.0,1.0),(3.0,3.0),(1.0,1.0)],
                               crs=LonLatWGS84)

        self.points3 = [Point((1, 1, 0), crs=LonLatWGS84),
                        Point((3, 1, 3), crs=LonLatWGS84),
                        Point((4, 3, 2), crs=LonLatWGS84),
                        Point((2, 2, -1), crs=LonLatWGS84)]

        self.line3 = Line([(1,5,2),(5,5,-1),(5,1,3),(3,3,1),(1,1,0)], crs=LonLatWGS84)

        self.polygon3 = Polygon([(1,5,2),(5,5,-1),(5,1,3),(3,3,1),(1,1,0)], crs=LonLatWGS84)

        testfiles = ["points.shp", "line.shp", "polygon.shp"]
        if any(not exists(join(TMPDATA, "shapefiles/", fnm)) for fnm in testfiles):
            self.saveTestData()
        return

    def saveTestData(self):
        testfiles = [(self.multipoint, "points"),
                     (self.line, "line"),
                     (self.polygon, "polygon")]
        os.makedirs(os.path.join(TMPDATA, "shapefiles"))
        for (geom, fnm) in testfiles:
            geom.to_shapefile(os.path.join(TMPDATA, "shapefiles", fnm))
        return

    def assertGeomEqual(self, this, that):
        self.assertTrue(np.all(this.get_vertices() == that.get_vertices()))
        try:
            self.assertEqual(this.crs.get_proj4(), that.crs.get_proj4())
        except AttributeError:
            print("warning: crs equality not established")
        return

    def test_writepoint(self):
        point = self.points[0]
        point.to_shapefile(os.path.join(TESTDIR, "data/point"))
        for fnm in ("point.shx", "point.shx", "point.dbf", "point.prj"):
            self.assertTrue(os.path.isfile(os.path.join(TESTDIR, "data", fnm)))
        return

    def test_writepoints(self):
        points = self.points
        shp.write_shapefile(os.path.join(TESTDIR, "data/points.shp"), *points)
        for fnm in ("points.shx", "points.shx", "points.dbf", "points.prj"):
            self.assertTrue(os.path.isfile(os.path.join(TESTDIR, "data", fnm)))
        return

    def test_writemultipoint(self):
        mp = Multipoint(self.points)
        mp.to_shapefile(os.path.join(TESTDIR, "data/multipoint"))
        for fnm in ("multipoint.shx", "multipoint.shx", "multipoint.dbf", "multipoint.prj"):
            self.assertTrue(os.path.isfile(os.path.join(TESTDIR, "data", fnm)))
        return

    def test_writeline(self):
        self.line.to_shapefile(os.path.join(TESTDIR, "data/line"))
        for fnm in ("line.shx", "line.shx", "line.dbf", "line.prj"):
            self.assertTrue(os.path.isfile(os.path.join(TESTDIR, "data", fnm)))
        return

    def test_writepoly(self):
        self.polygon.to_shapefile(os.path.join(TESTDIR, "data/polygon"))
        for fnm in ("polygon.shx", "polygon.shx", "polygon.dbf", "polygon.prj"):
            self.assertTrue(os.path.isfile(os.path.join(TESTDIR, "data", fnm)))
        return

    def test_writepoints3(self):
        mp = Multipoint(self.points3)
        mp.to_shapefile(os.path.join(TESTDIR, "data/multipointz"))
        for fnm in ("multipointz.shx", "multipointz.shx", "multipointz.dbf", "multipointz.prj"):
            self.assertTrue(os.path.isfile(os.path.join(TESTDIR, "data", fnm)))
        return

    def test_writeline3(self):
        self.line3.to_shapefile(os.path.join(TESTDIR, "data/linez"))
        for fnm in ("linez.shx", "linez.shx", "linez.dbf", "linez.prj"):
            self.assertTrue(os.path.isfile(os.path.join(TESTDIR, "data", fnm)))
        return

    def test_writepoly3(self):
        self.polygon3.to_shapefile(os.path.join(TESTDIR, "data/polygonz"))
        for fnm in ("polygonz.shx", "polygonz.shx", "polygonz.dbf", "polygonz.prj"):
            self.assertTrue(os.path.isfile(os.path.join(TESTDIR, "data", fnm)))
        return

    def test_write_collection_multipoint(self):
        mp = Multipoint([p.vertex for p in self.points])
        mp0 = copy(mp)
        mp1 = copy(mp.shift((4, 2)))
        mp2 = copy(mp.shift((-2, 3)))
        shp.write_shapefile(os.path.join(TESTDIR, "data/mp_collection.shp"),
                            mp0, mp1, mp2)
        for fnm in ("mp_collection.shx", "mp_collection.shx", "mp_collection.dbf", "mp_collection.prj"):
            self.assertTrue(os.path.isfile(os.path.join(TESTDIR, "data", fnm)))
        return

    def test_write_collection_lines(self):
        line0 = copy(self.line)
        line1 = copy(self.line.shift((4, 2)))
        line2 = copy(self.line.shift((-2, 3)))
        shp.write_shapefile(os.path.join(TESTDIR, "data/line_collection.shp"),
                            line0, line1, line2)
        for fnm in ("line_collection.shx", "line_collection.shx", "line_collection.dbf", "line_collection.prj"):
            self.assertTrue(os.path.isfile(os.path.join(TESTDIR, "data", fnm)))
        return

    def test_read_points(self):
        points = read_shapefile(os.path.join(TESTDATA, "shp_input", "points"))
        self.assertEqual(len(points), 4)
        pt = points[0]
        self.assertTrue("+proj=lonlat" in pt.crs.get_proj4())
        self.assertTrue("+a=6378137.0" in pt.crs.get_proj4())
        self.assertTrue("+f=0.00335281" in pt.crs.get_proj4())
        mp = Multipoint(points)
        self.assertEqual(mp.d["species"], ['T. officianale', 'C. tectorum', 'M. alba', 'V. cracca'])
        self.assertEqual(mp.d["ID"], ['0', '1', '2', '3'])
        self.assertEqual(mp.coordinates, ((1.0, 3.0, 4.0, 2.0), (1.0, 1.0, 3.0, 2.0)))

    def test_read_line(self):
        line = read_shapefile(os.path.join(TESTDATA, "shp_input", "line"))[0]
        self.assertTrue("+proj=lonlat" in line.crs.get_proj4())
        self.assertTrue("+a=6378137.0" in line.crs.get_proj4())
        self.assertTrue("+f=0.00335281" in line.crs.get_proj4())
        self.assertEqual(line.coordinates, ((1.0, 5.0, 5.0, 3.0, 1.0), (5.0, 5.0, 1.0, 3.0, 1.0)))
        return

    def test_read_polygon(self):
        polygon = read_shapefile(os.path.join(TESTDATA, "shp_input", "polygon"))[0]
        self.assertTrue("+proj=lonlat" in polygon.crs.get_proj4())
        self.assertTrue("+a=6378137.0" in polygon.crs.get_proj4())
        self.assertTrue("+f=0.00335281" in polygon.crs.get_proj4())
        self.assertEqual(polygon.coordinates, ((1.0, 5.0, 5.0, 3.0, 1.0), (5.0, 5.0, 1.0, 3.0, 1.0)))
        return

    def test_read_points_newp(self):
        # Read a multipoint with a projected cooridnate system
        newp = read_shapefile(os.path.join(TESTDATA, "shp_input", 
                                                     "newp_nsidc_north"))

        proj4 = ('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 '
                 '+y_0=0 +a=6378273 +b=6356889.449 +units=m +no_defs')

        for part in proj4.split():
            self.assertTrue(part[:8] in newp[0].crs.get_proj4())

        coords = list(zip(*[pt.vertex[:2] for pt in newp]))
        self.assertEqual(coords, [(521236.8297444395, 521236.8297444395,
                                   521236.8297444395, 547490.4452879033,
                                   547490.4452879033, 547490.4452879033,
                                   587584.1578033275, 587584.1578033275,
                                   587584.1578033275, 571828.4918982167,
                                   571828.4918982167),
                                  (-888853.1384770898, -888853.1384770898,
                                   -888853.1384770898, -902049.3617542256,
                                   -902049.3617542256, -902049.3617542256,
                                   -871214.0673764511, -871214.0673764511,
                                   -871214.0673764511, -850080.914674058,
                                   -850080.914674058)])

        meterno = [pt.properties["meterno"] for pt in newp]
        self.assertEqual(meterno, ['IMS1/1', 'IMS2/1', '5952/2', 'IMS4/1',
                                   '5953/2', '1963/13', 'IMS5/1', '5213/A',
                                   '2121/13', 'IMS3/1', '3613/2'])

        depth = [pt.properties["depth_m"] for pt in newp]
        self.assertEqual(depth, ['73', '143', '247', '86', '147', '250', '74', 
                                 '142', '235', '150', '248'])
        return

class ShapefileAttributeTests(unittest.TestCase):

    def test_infer_ogr_fieldtype(self):
        self.assertEqual(shp.ogr_get_fieldtype(1), (0, 32))
        self.assertEqual(shp.ogr_get_fieldtype([1, 2]), (1, 1000))

        self.assertEqual(shp.ogr_get_fieldtype(1.0), (2, 32))
        self.assertEqual(shp.ogr_get_fieldtype([1.0, 1.5]), (3, 1000))

        # everything should be interpretted as WideString
        self.assertEqual(shp.ogr_get_fieldtype("hello"), (4, 180))
        self.assertEqual(shp.ogr_get_fieldtype(["list","of","strings"]),(5, 1000))

        # doesn't work on Python 2
        #self.assertEqual(shp.ogr_get_fieldtype(b'0b110001'), 8)

        # dates
        self.assertEqual(shp.ogr_get_fieldtype(datetime.date(2013, 11, 17)), (9, 32))
        self.assertEqual(shp.ogr_get_fieldtype(datetime.time(8, 30, 0)), (10, 32))
        self.assertEqual(shp.ogr_get_fieldtype(datetime.datetime(2013, 11, 17, 8, 30, 0)), (11, 64))
        return


class ShapelibTestSuite(unittest.TestCase):
    """ Open and verify the shapefiles provided with the shapelib testsuite. """

    def setUp(self):
        self.dirname = os.path.join(TESTDATA, "shapefile", "shapelib")

    def test_(self):
        res = read_shapefile(os.path.join(self.dirname, "test.shp"))

    def test_0(self):
        res = read_shapefile(os.path.join(self.dirname, "test0.shp"))

    def test_1(self):
        res = read_shapefile(os.path.join(self.dirname, "test1.shp"))

    def test_2(self):
        res = read_shapefile(os.path.join(self.dirname, "test2.shp"))

    def test_3(self):
        res = read_shapefile(os.path.join(self.dirname, "test3.shp"))

    def test_4(self):
        res = read_shapefile(os.path.join(self.dirname, "test4.shp"))

    def test_5(self):
        res = read_shapefile(os.path.join(self.dirname, "test5.shp"))

    def test_6(self):
        res = read_shapefile(os.path.join(self.dirname, "test6.shp"))

    def test_7(self):
        res = read_shapefile(os.path.join(self.dirname, "test7.shp"))

    def test_8(self):
        res = read_shapefile(os.path.join(self.dirname, "test8.shp"))

    def test_9(self):
        res = read_shapefile(os.path.join(self.dirname, "test9.shp"))

    def test_10(self):
        res = read_shapefile(os.path.join(self.dirname, "test10.shp"))

    def test_11(self):
        res = read_shapefile(os.path.join(self.dirname, "test11.shp"))

    def test_12(self):
        res = read_shapefile(os.path.join(self.dirname, "test12.shp"))

    def test_13(self):
        res = read_shapefile(os.path.join(self.dirname, "test13.shp"))

if __name__ == "__main__":
    unittest.main()
