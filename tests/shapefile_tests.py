import unittest
import os
from os.path import exists, join
import numpy as np
from test_helper import TESTDIR, TESTDATA
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
        if any(not exists(join(TESTDATA, "shapefiles/", fnm)) for fnm in testfiles):
            self.saveTestData()
        return

    def saveTestData(self):
        testfiles = [(self.multipoint, "points"),
                     (self.line, "line"),
                     (self.polygon, "polygon")]
        os.makedirs(os.path.join(TESTDATA, "shapefiles"))
        for (geom, fnm) in testfiles:
            geom.to_shapefile(os.path.join(TESTDATA, "shapefiles", fnm))
        return

    def assertGeomEqual(self, this, that):
        self.assertTrue(np.all(this.get_vertices() == that.get_vertices()))
        try:
            self.assertEqual(this.crs.get_proj4(), that.crs.get_proj4())
        except AttributeError:
            print("warning: crs equality not established")
        return

    def test_writepoints(self):
        mp = Multipoint(self.points)
        mp.to_shapefile(os.path.join(TESTDIR, "data/points_shp"))
        for fnm in ("points_shp.shx", "points_shp.shx", "points_shp.dbf", "points_shp.prj"):
            self.assertTrue(os.path.isfile(os.path.join(TESTDIR, "data", fnm)))
        return

    def test_writeline(self):
        self.line.to_shapefile(os.path.join(TESTDIR, "data/line_shp"))
        for fnm in ("line_shp.shx", "line_shp.shx", "line_shp.dbf", "line_shp.prj"):
            self.assertTrue(os.path.isfile(os.path.join(TESTDIR, "data", fnm)))
        return

    def test_writepoly(self):
        self.polygon.to_shapefile(os.path.join(TESTDIR, "data/polygon_shp"))
        for fnm in ("polygon_shp.shx", "polygon_shp.shx", "polygon_shp.dbf", "polygon_shp.prj"):
            self.assertTrue(os.path.isfile(os.path.join(TESTDIR, "data", fnm)))
        return

    def test_writepoints3(self):
        mp = Multipoint(self.points3)
        mp.to_shapefile(os.path.join(TESTDIR, "data/pointsz_shp"))
        for fnm in ("pointsz_shp.shx", "pointsz_shp.shx", "pointsz_shp.dbf", "pointsz_shp.prj"):
            self.assertTrue(os.path.isfile(os.path.join(TESTDIR, "data", fnm)))
        return

    def test_writeline3(self):
        self.line3.to_shapefile(os.path.join(TESTDIR, "data/linez_shp"))
        for fnm in ("linez_shp.shx", "linez_shp.shx", "linez_shp.dbf", "linez_shp.prj"):
            self.assertTrue(os.path.isfile(os.path.join(TESTDIR, "data", fnm)))
        return

    def test_writepoly3(self):
        self.polygon3.to_shapefile(os.path.join(TESTDIR, "data/polygonz_shp"))
        for fnm in ("polygonz_shp.shx", "polygonz_shp.shx", "polygonz_shp.dbf", "polygonz_shp.prj"):
            self.assertTrue(os.path.isfile(os.path.join(TESTDIR, "data", fnm)))
        return

    def test_write_collection_points(self):
        mp = Multipoint([p.vertex for p in self.points])
        mp0 = copy(mp)
        mp1 = copy(mp.shift((4, 2)))
        mp2 = copy(mp.shift((-2, 3)))
        shp.write_shapefile(os.path.join(TESTDIR, "data/points_collection.shp"),
                            mp0, mp1, mp2)
        for fnm in ("points_collection.shx", "points_collection.shx", "points_collection.dbf", "points_collection.prj"):
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

if __name__ == "__main__":
    unittest.main()
