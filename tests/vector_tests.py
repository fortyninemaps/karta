""" Unit tests for vector functions """

import unittest
import os
import math
import numpy as np
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import json
from test_helper import md5sum, md5sum_file, TESTDATA, TESTDIR

import karta.vector as vector
import karta.crs as crs
from karta.vector.geojson import GeoJSONReader
from karta.vector.guppy import Point, Multipoint, Line, Polygon
from karta.vector.metadata import Metadata

class TestGuppy(unittest.TestCase):

    def setUp(self):
        self.point = Point((1.0, 2.0, 3.0), data={"color":(43,67,10)},
                           properties="apple")

        self.vertices = [(2.0, 9.0, 9.0), (4.0, 1.0, 9.0), (4.0, 1.0, 5.0),
                         (2.0, 8.0, 0.0), (9.0, 8.0, 4.0), (1.0, 4.0, 6.0),
                         (7.0, 3.0, 4.0), (2.0, 5.0, 3.0), (1.0, 6.0, 6.0),
                         (8.0, 1.0, 0.0), (5.0, 5.0, 1.0), (4.0, 5.0, 7.0),
                         (3.0, 3.0, 5.0), (9.0, 0.0, 9.0), (6.0, 3.0, 8.0),
                         (4.0, 5.0, 7.0), (9.0, 9.0, 4.0), (1.0, 4.0, 7.0),
                         (1.0, 7.0, 8.0), (9.0, 1.0, 6.0)]

        self.data = [99.0, 2.0, 60.0, 75.0, 71.0, 34.0, 1.0, 49.0, 4.0, 36.0,
                     47.0, 58.0, 65.0, 72.0, 4.0, 27.0, 52.0, 37.0, 95.0, 17.0]

        self.mp = Multipoint(self.vertices, data=self.data)
        self.line = Line(self.vertices, data=self.data)
        self.poly = Polygon([(0.0, 8.0), (0.0, 5.0), (6.0, 1.0)])
        self.ring = Polygon([(2.0, 2.0), (4.0, 2.0), (3.0, 6.0)])
        self.ringed_poly = Polygon([(0.0, 0.0), (10, 0.0),
                                    (10.0, 10.0), (0.0, 10.0)],
                                   subs=[self.ring])
        self.unitsquare = Polygon([(0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,1.0)])
        return


    def test_point_equality(self):
        pt1 = Point((3.0, 4.0))
        pt2 = Point((3.0, 4.0, 5.0))
        pt3 = Point((3.0, 4.0, 5.0), data={"species":"T. officianale", "density":"high"})
        self.assertFalse(pt1 == pt2)
        self.assertFalse(pt1 == pt3)
        self.assertFalse(pt2 == pt3)
        return

    def test_point_vertex(self):
        self.assertEqual(self.point.get_vertex(), (1.0, 2.0, 3.0))
        return

    def test_point_coordsxy(self):
        self.assertEqual(self.point.coordsxy(), (1.0, 2.0))
        self.assertEqual(self.point[0], 1.0)
        self.assertEqual(self.point[1], 2.0)
        return

    def test_point_azimuth(self):
        point = Point((1.0, 2.0))

        other = Point((2.0, 3.0))
        self.assertEqual(point.azimuth(other), 0.25*np.pi)

        other = Point((0.0, 3.0))
        self.assertEqual(point.azimuth(other), 1.75*np.pi)

        other = Point((0.0, 1.0))
        self.assertEqual(point.azimuth(other), 1.25 * np.pi)

        other = Point((2.0, 1.0))
        self.assertEqual(point.azimuth(other), 0.75 * np.pi)

        other = Point((1.0, 3.0))
        self.assertEqual(point.azimuth(other), 0.0)

        other = Point((1.0, 1.0))
        self.assertEqual(point.azimuth(other), np.pi)
        return

    def test_point_azimuth2(self):
        point = Point((5.0, 2.0))
        other = Point((5.0, 2.0))
        self.assertTrue(np.isnan(point.azimuth(other)))
        return

    def test_point_shift(self):
        point = Point((-3.0, 5.0, 2.5), data={"color":(43,67,10)},
                      properties="apple")
        point.shift((4.0, -3.0, 0.5))
        self.assertEqual(self.point, point)
        return

    def test_nearest_to(self):
        self.assertEqual(self.mp.nearest_point_to(self.point), self.mp[12])
        return

    def test_multipoint_zip_init(self):
        x = range(-10, 10)
        y = [_x**2 for _x in x]
        Line(zip(x, y))
        return

    def test_multipoint_shift(self):
        vertices = [(a-1,b+2,c-0.5) for (a,b,c) in self.vertices]
        mp = Multipoint(vertices, data=self.data)
        mp.shift((1, -2, 0.5))
        self.assertEqual(mp, self.mp)
        return

    def test_multipoint_subset(self):
        ss1 = self.mp._subset(range(2,7))
        ss2 = self.line._subset(range(2,7))
        self.assertTrue(isinstance(ss1, Multipoint))
        self.assertTrue(isinstance(ss2, Line))
        return

    def test_multipoint_get(self):
        self.assertEqual(self.mp[0], Point(self.vertices[0],
                                           data=self.mp.data[0],
                                           properties=self.mp.properties))
        return

    def test_multipoint_set(self):
        mp1 = Multipoint([(3.0, 3.0), (5.0, 1.0), (3.0, 1.0),
                         (4.0, 4.0), (0.0, 1.0)],
                         data=["rankin", "corbet", "arviat",
                               "severn", "churchill"])
        mp2 = Multipoint([(3.0, 3.0), (5.0, 1.0), (4.0, 5.0),
                         (4.0, 4.0), (0.0, 1.0)],
                         data=["rankin", "corbet", "umiujaq",
                               "severn", "churchill"])
        mp1[2] = (4.0, 5.0)
        self.assertNotEqual(mp1, mp2)
        mp1[2] = Point((4.0, 5.0), data="umiujaq")
        self.assertEqual(mp1, mp2)
        return

    def test_multipoint_slicing(self):
        submp = Multipoint(self.vertices[5:10], data=self.data[5:10])
        self.assertEqual(self.mp[5:10], submp)

        submp = Multipoint(self.vertices[5:], data=self.data[5:])
        self.assertEqual(self.mp[5:], submp)
        return

    def test_multipoint_negative_index(self):
        self.assertEqual(self.mp[len(self.mp)-1], self.mp[-1])
        return

    def test_multipoint_bbox(self):
        bbox = (1.0, 9.0, 0.0, 9.0, 0.0, 9.0)
        self.assertEqual(self.mp.get_bbox(), bbox)
        return

    def test_multipoint_bbox_overlap(self):
        self.assertTrue(self.mp._bbox_overlap(self.poly))
        return

    def test_multipoint_datadict(self):
        # create a line
        vertices = [(2.0, 9.0, 9.0), (4.0, 1.0, 9.0), (4.0, 1.0, 5.0),
                    (2.0, 8.0, 0.0), (9.0, 8.0, 4.0), (1.0, 4.0, 6.0),
                    (7.0, 3.0, 4.0), (2.0, 5.0, 3.0), (1.0, 6.0, 6.0),
                    (8.0, 1.0, 0.0), (5.0, 5.0, 1.0), (4.0, 5.0, 7.0),
                    (3.0, 3.0, 5.0), (9.0, 0.0, 9.0), (6.0, 3.0, 8.0),
                    (4.0, 5.0, 7.0), (9.0, 9.0, 4.0), (1.0, 4.0, 7.0),
                    (1.0, 7.0, 8.0), (9.0, 1.0, 6.0)]

        data0 = [99.0, 2.0, 60.0, 75.0, 71.0, 34.0, 1.0, 49.0, 4.0, 36.0, 47.0,
                 58.0, 65.0, 72.0, 4.0, 27.0, 52.0, 37.0, 95.0, 17.0]

        data1 = [54.0, 40.0, 77.0, 18.0, 84.0, 91.0, 61.0, 92.0, 19.0, 42.0,
                 50.0, 25.0, 11.0, 80.0, 59.0, 56.0, 32.0, 8.0, 88.0, 76.0]

        L2 = Multipoint(vertices, data={'d0':data0, 'd1':data1})
        return

    def test_multipoint_within_radius(self):
        vertices = [(float(x),float(y)) for x in range(-10,11)
                                        for y in range(-10,11)]
        ans = [v for v in vertices if math.sqrt(v[0]**2 + v[1]**2) <= 5.0]
        mp = Multipoint(vertices)
        sub = mp.within_radius(Point((0,0)), 5.0)
        self.assertEqual(sub, Multipoint(ans))
        return

    def test_multipoint_within_bbox(self):
        vertices = [(float(x),float(y)) for x in range(-10,11)
                                        for y in range(-10,11)]
        ans = [v for v in vertices if (-5.0<=v[0]<=5.0) and (-4.0<=v[1]<=6.0)]
        mp = Multipoint(vertices)
        sub = mp.within_bbox((-5.0, 5.0, -4.0, 6.0))
        self.assertEqual(sub, Multipoint(ans))
        return

    def test_multipoint_from_points(self):
        pts = [Point((x,y), data={"amount": x+y}) for x in range(3)
                                                  for y in range(2)]
        mp = Multipoint(pts)
        verts = [Point((x,y), data={"amount": x+y}) for x in range(3)
                                                    for y in range(2)]
        ans = Multipoint(verts, data={"amount":[sum(v) for v in verts]})
        self.assertEqual(mp, ans)
        return

    def test_connected_multipoint_shortest_distance_to(self):
        line = Line([(0.0, 0.0), (2.0, 2.0), (5.0, 4.0)])
        dist = line.shortest_distance_to(Point((0.0, 2.0)))
        self.assertTrue(abs(dist - math.sqrt(2)) < 1e-10)
        return

    def test_connected_multipoint_shortest_distance_to2(self):
        line = Line([(127.0, -35.0), (132.0, -28.0), (142.0, -29.0)], crs=crs.LONLAT)
        dist = line.shortest_distance_to(Point((98.0, -7.0), crs=crs.LONLAT))
        self.assertAlmostEqual(dist, 4257313.5324397, places=6)
        return

    def test_connected_multipoint_nearest_on_boundary(self):
        line = Line([(0.0, 0.0), (2.0, 2.0), (5.0, 4.0)])
        npt = line.nearest_on_boundary(Point((0.0, 2.0)))
        self.assertEqual(npt, Point((1.0, 1.0)))
        return

    def assertPointAlmostEqual(self, a, b):
        for (a_, b_) in zip(a.vertex, b.vertex):
            self.assertAlmostEqual(a_, b_, places=5)
        self.assertEqual(a.data, b.data)
        self.assertEqual(a.properties, b.properties)
        self.assertEqual(a._crs, b._crs)
        return

    def test_connected_multipoint_nearest_on_boundary2(self):
        line = Line([(-40, 0.0), (35, 0.0)], crs=crs.LONLAT)
        npt = line.nearest_on_boundary(Point((30.0, 80.0), crs=crs.LONLAT))
        self.assertPointAlmostEqual(npt, Point((30.0, 0.0), crs=crs.LONLAT))
        return

    # def test_connected_multipoint_nearest_on_boundary3(self):
    #     # This is the test that tends to break naive root finding schemes
    #     line = Line([(-40, 0.0), (35, 0.0)], crs=crs.LONLAT)
    #     npt = line.nearest_on_boundary(Point((30.0, 1e-8), crs=crs.LONLAT))
    #     self.assertPointAlmostEqual(npt, Point((30.0, 0.0), crs=crs.LONLAT))
    #     return

    # def test_connected_multipoint_nearest_on_boundary4(self):
    #     line = Line([(-20.0, 32.0), (-26.0, 43.0), (-38.0, 39.0)], crs=crs.LONLAT)
    #     npt = line.nearest_on_boundary(Point((-34.0, 52.0), crs=crs.LONLAT))
    #     self.assertPointAlmostEqual(npt, Point((-27.98347, 42.456316), crs=crs.LONLAT))
    #     return

    def test_line_add_vertex2d(self):
        ln0 = Line([(3.0, 3.0), (5.0, 1.0), (3.0, 1.0)])
        ln1 = Line([(3.0, 3.0), (5.0, 1.0), (3.0, 1.0),
                    (4.0, 4.0), (0.0, 1.0)])
        ln0.add_vertex((4.0, 4.0))
        ln0.add_vertex(Point((0.0, 1.0)))
        self.assertEqual(ln0, ln1)
        return

    def test_line_add_vertex3d(self):
        ln0 = Line([(3.0, 3.0, 2.0), (5.0, 1.0, 0.0), (3.0, 1.0, 5.0)])
        ln1 = Line([(3.0, 3.0, 2.0), (5.0, 1.0, 0.0), (3.0, 1.0, 5.0),
                    (4.0, 4.0, 6.0), (0.0, 1.0, 3.0)])
        ln0.add_vertex((4.0, 4.0, 6.0))
        ln0.add_vertex(Point((0.0, 1.0, 3.0)))
        self.assertEqual(ln0, ln1)
        return

    def test_line_extend(self):
        ln0a = Line([(3.0, 3.0, 2.0), (5.0, 1.0, 0.0), (3.0, 1.0, 5.0)])
        ln0b = Line([(4.0, 4.0, 6.0), (0.0, 1.0, 3.0)])
        ln1 = Line([(3.0, 3.0, 2.0), (5.0, 1.0, 0.0), (3.0, 1.0, 5.0),
                    (4.0, 4.0, 6.0), (0.0, 1.0, 3.0)])
        ln0a.extend(ln0b)
        self.assertEqual(ln0a, ln1)

    def test_line_remove(self):
        ln = Line([(3.0, 3.0, 2.0), (5.0, 1.0, 0.0), (3.0, 1.0, 5.0),
                   (4.0, 4.0, 6.0), (0.0, 1.0, 3.0)],
                  data=["red", "green", "blue", "chartreuse", "aquamarine"])
        lnresult = Line([(3.0, 3.0, 2.0), (5.0, 1.0, 0.0), (3.0, 1.0, 5.0),
                         (0.0, 1.0, 3.0)],
                        data=["red", "green", "blue", "aquamarine"])
        pt = ln.remove_vertex(3)
        ptresult = Point((4.0, 4.0, 6.0), data="chartreuse")
        self.assertEqual(pt, ptresult)
        self.assertEqual(ln, lnresult)
        return

    def test_line_intersection(self):
        line0 = Line([(0.0, 0.0), (3.0, 3.0)])
        line1 = Line([(0.0, 3.0), (3.0, 0.0)])
        self.assertTrue(line0.intersects(line1))
        self.assertEqual(line0.intersections(line1), Multipoint([(1.5, 1.5)]))
        return

    def test_poly_vertices(self):
        self.assertTrue((self.poly.get_vertices() ==
            np.array([(0.0, 8.0), (0.0, 5.0), (6.0, 1.0), (0.0, 8.0)])).all())
        return

    def test_poly_coordinates(self):
        self.assertEqual(self.poly.get_coordinate_lists(),
                         ([0.0, 0.0, 6.0, 0.0], [8.0, 5.0, 1.0, 8.0]))
        return

    def test_poly_extents(self):
        self.assertEqual(self.poly.get_extents(), [(0.0, 6.0), (1.0, 8.0)])
        return

    def test_poly_length(self):
        self.assertEqual(self.poly.length(), 19.430647008220866)
        return

    def test_poly_contains(self):
        pt0 = Point((-0.5, 0.92))
        pt1 = Point((0.125, 0.875))
        self.assertFalse(self.unitsquare.contains(pt0))
        self.assertTrue(self.unitsquare.contains(pt1))
        return

    def test_poly_getitem(self):
        poly = Polygon([(0.0, 8.0), (0.0, 5.0), (6.0, 1.0), (7.0, 2.0),
                        (5.0, 4.0)])
        sub = poly[:3]
        self.assertEqual(sub, Line([(0.0, 8.0), (0.0, 5.0), (6.0, 1.0)]))
        return

    def test_poly_getitem2(self):
        poly = Polygon([(0.0, 8.0), (0.0, 5.0), (6.0, 1.0), (7.0, 2.0),
                        (5.0, 4.0)])
        sub = poly[:4:2]
        self.assertEqual(sub, Line([(0.0, 8.0), (6.0, 1.0)]))
        return

    def test_poly_getitem3(self):
        poly = Polygon([(0.0, 8.0), (0.0, 5.0), (6.0, 1.0), (7.0, 2.0),
                        (5.0, 4.0)])
        sub = poly[:]
        self.assertEqual(sub, poly)
        return

    def test_ringedpoly_perimeter(self):
        self.assertEqual(round(self.ringed_poly.perimeter(), 3), 50.246)
        return

    def test_ringedpoly_area(self):
        self.assertEqual(self.ringed_poly.area(), 100 - self.ring.area())
        return

    def test_segments(self):
        v = self.vertices
        self.assertEqual([tuple(a.vertices) for a in self.line.segments()],
                         [(v[i], v[i+1]) for i in range(len(self.vertices)-1)])
        return

    def test_within_distance(self):
        line = Line([(0,0), (1,1), (3,1)])
        pt = Point((1,1.5))
        self.assertTrue(line.within_distance(pt, 0.6))
        self.assertFalse(line.within_distance(pt, 0.4))
        return

    def test_walk_cartesian(self):
        start = Point((-3, -4), crs=crs.CARTESIAN)
        dest = start.walk(5.0, math.atan(4.0/3.0))
        self.assertAlmostEqual(dest.x, 0.0)
        self.assertAlmostEqual(dest.y, 0.0)
        return


class TestGuppyProj(unittest.TestCase):

    def setUp(self):
        self.vancouver = Point((-123.1, 49.25), crs=crs.LONLAT)
        self.ottawa = Point((-75.69, 45.42), crs=crs.LONLAT)
        self.whitehorse = Point((-135.05, 60.72), crs=crs.LONLAT)
        return

    def test_greatcircle(self):
        d1 = self.vancouver.greatcircle(self.ottawa)
        d2 = self.vancouver.greatcircle(self.whitehorse)
        d3 = self.whitehorse.greatcircle(self.ottawa)
        d4 = self.whitehorse.greatcircle(self.vancouver)
        self.assertTrue(abs(d1 - 3549030.70541) < 1e-5)
        self.assertTrue(abs(d2 - 1483327.53922) < 1e-5)
        self.assertTrue(abs(d3 - 4151366.88185) < 1e-5)
        self.assertTrue(abs(d4 - 1483327.53922) < 1e-5)
        return

    def test_azimuth_lonlat(self):
        az1 = self.vancouver.azimuth(self.ottawa) * 180. / math.pi
        az2 = self.vancouver.azimuth(self.whitehorse) * 180. / math.pi
        self.assertTrue(abs(az1 - 78.4833443403) < 1e-10)
        self.assertTrue(abs(az2 - (-26.1358265354)) < 1e-10)
        return

    def test_walk_lonlat(self):
        start = Point((-132.14, 54.01), crs=crs.LONLAT_WGS84)
        dest = start.walk(5440.0, 106.8*math.pi/180.0)
        self.assertAlmostEqual(dest.x, -132.0605910876)
        self.assertAlmostEqual(dest.y, 53.99584742821)
        return


class TestGuppyOutput(unittest.TestCase):

    def setUp(self):

        vertices = [(2.0, 9.0, 9.0), (4.0, 1.0, 9.0), (4.0, 1.0, 5.0),
                    (2.0, 8.0, 0.0), (9.0, 8.0, 4.0), (1.0, 4.0, 6.0),
                    (7.0, 3.0, 4.0), (2.0, 5.0, 3.0), (1.0, 6.0, 6.0),
                    (8.0, 1.0, 0.0), (5.0, 5.0, 1.0), (4.0, 5.0, 7.0),
                    (3.0, 3.0, 5.0), (9.0, 0.0, 9.0), (6.0, 3.0, 8.0),
                    (4.0, 5.0, 7.0), (9.0, 9.0, 4.0), (1.0, 4.0, 7.0),
                    (1.0, 7.0, 8.0), (9.0, 1.0, 6.0)]

        data0 = [99.0, 2.0, 60.0, 75.0, 71.0, 34.0, 1.0, 49.0, 4.0, 36.0, 47.0,
                 58.0, 65.0, 72.0, 4.0, 27.0, 52.0, 37.0, 95.0, 17.0]

        data1 = [54.0, 40.0, 77.0, 18.0, 84.0, 91.0, 61.0, 92.0, 19.0, 42.0,
                 50.0, 25.0, 11.0, 80.0, 59.0, 56.0, 32.0, 8.0, 88.0, 76.0]
        self.mp = Multipoint(vertices, data={'d0':data0, 'd1':data1})

    # Due to the output of ElementTree.tostring not being deterministic, this
    # test randomly fails due to the DataArray attributes being swapped. The
    # output is correct, but it doesn't match the reference data. This is a bug
    # in the test.
    #def test_mp2vtp(self):
    #    # Test VTK output for a Multipoint
    #    s = StringIO()
    #    self.mp.to_vtk(s)
    #    s.seek(0)
    #    #a1 = md5sum(s)
    #    #a2 = md5sum_file(os.path.join(TESTDATA, 'testmp2vtp.vtp'))
    #    #print(a1)
    #    #print(a2)
    #    #self.assertEqual(md5sum(s),
    #    #                 md5sum_file(os.path.join(TESTDATA, 'testmp2vtp.vtp')))
    #    return

class TestMetadata(unittest.TestCase):

    def setUp(self):
        self.onefield = Metadata(data=np.arange(200), singleton=False)
        self.multifield = Metadata(
                            data={"a":np.arange(200), "b":np.arange(200,400), 
                                  "c":np.arange(200)**2},
                            singleton=False)
        return

    def test_indexing(self):
        self.assertEqual(self.onefield[10], {"values": 10})
        self.assertEqual(self.multifield[10], {"a":10, "b": 210, "c": 100})
        return

    def test_slicing(self):
        self.assertTrue(np.all(self.onefield[5:15]["values"] == np.arange(5, 15)))
        res = self.multifield[5:15]
        self.assertTrue(np.all(res["a"] == np.arange(5, 15)))
        self.assertTrue(np.all(res["b"] == np.arange(205, 215)))
        self.assertTrue(np.all(res["c"] == np.arange(5, 15)**2))
        return

    def test_keys(self):
        self.assertTrue(np.all(self.onefield["values"] == np.arange(200)))
        self.assertTrue(np.all(self.multifield["a"] == np.arange(200)))
        self.assertTrue(np.all(self.multifield["b"] == np.arange(200, 400)))
        self.assertTrue(np.all(self.multifield["c"] == np.arange(200)**2))
        return

    def test_reference_vs_value_list(self):
        # metadata objects should carry values rather than references
        L = [1,2,3,4,5]
        md = Metadata(L)
        L[3] = -99
        self.assertEqual(md["values"], [1,2,3,4,5])
        return

    def test_reference_vs_value_dict(self):
        # metadata objects should carry values rather than references
        D = {"A":[1,2,3,4,5],
             "B":[6,7,8,9,10]}
        md = Metadata(D)
        D["A"][3] = -99
        self.assertEqual(md["A"], [1,2,3,4,5])
        return

if __name__ == "__main__":
    unittest.main()

