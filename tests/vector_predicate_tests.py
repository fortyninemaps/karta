""" Unit tests for vector geometry predicate methods """

from __future__ import division
import unittest
import numpy as np

from karta.vector.geometry import (Point, Line, Polygon,
                                   Multipoint, Multiline, Multipolygon)
from karta.crs import (Cartesian, SphericalEarth, LonLatWGS84)
from karta.errors import CRSError

class TestUnaryPredicates(unittest.TestCase):

    def test_poly_clockwise(self):
        p = Polygon([(0,0), (0,1), (1,1), (1,0)])
        self.assertTrue(p.isclockwise())
        return

    def test_poly_counterclockwise(self):
        p = Polygon([(0,0), (1,0), (1,1), (0,1)])
        self.assertFalse(p.isclockwise())
        return

    def test_poly_polar(self):
        p = Polygon([(0.0, 80.0), (30.0, 80.0), (60.0, 80.0), (90.0, 80.0),
                     (120.0, 80.0), (150.0, 80.0), (180.0, 80.0),
                     (-150.0, 80.0), (-120.0, 80.0), (-90.0, 80.0),
                     (-60.0, 80.0), (-30.0, 80.0)], crs=SphericalEarth)
        self.assertTrue(p.ispolar())

        p = Polygon([(0.0, 85.0, 0.0), (90.0, 85.0, 0.0), (180.0, 85.0, 0.0),
                     (-90.0, 85.0, 0.0)], crs=SphericalEarth)
        self.assertTrue(p.ispolar())

        p = Polygon([(45.0, 30.0), (40.0, 25.0), (45.0, 20.0), (35.0, 25.0)],
                    crs=SphericalEarth)
        self.assertFalse(p.ispolar())

        p = Polygon([(-80, 0), (-50, -10), (20, -8), (35, -17), (55, 15),
                     (-45, 18), (-60, 12)], crs=LonLatWGS84)
        self.assertFalse(p.ispolar())

        p = Polygon([(45.0, 30.0), (40.0, 25.0), (45.0, 20.0), (35.0, 25.0)],
                    crs=Cartesian)
        self.assertRaises(CRSError, p.ispolar)
        return

class TestBinaryPredicates(unittest.TestCase):

    def test_line_intersection(self):
        line0 = Line([(0.0, 0.0), (3.0, 3.0)])
        line1 = Line([(0.0, 3.0), (3.0, 0.0)])
        self.assertTrue(line0.intersects(line1))
        self.assertEqual(line0.intersections(line1), Multipoint([(1.5, 1.5)]))
        return

    def test_line_intersection2(self):
        # test lines that have overlapping bounding boxes, but don't cross
        #   -----
        #   | -----
        #   |     |
        #   ----- |
        #     -----
        line0 = Line([(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0)])
        line1 = Line([(1.0, 4.0), (-2.0, 4.0), (-2.0, 1.0), (1.0, 1.0)])
        self.assertFalse(line0.intersects(line1))
        return

    def test_poly_intersection(self):
        # test polygons formed exactly as in test_line_intersection2, except
        # the rings are implicitly closed
        #   -----
        #   | --x--
        #   | . . |
        #   --x-- |
        #     -----
        poly0 = Polygon([(0.0, 0.0), (3.0, 0.0), (3.0, 3.0), (0.0, 3.0)])
        poly1 = Polygon([(1.0, 4.0), (-2.0, 4.0), (-2.0, 1.0), (1.0, 1.0)])
        self.assertTrue(poly0.intersects(poly1))
        self.assertEqual(poly0.intersections(poly1), Multipoint([(0.0, 1.0), (1.0, 3.0)]))
        return

    def test_line_intersection_horizontal(self):
        line0 = Line([(-2.5, 2.5), (2.5, 2.5)])
        line1 = Line([(0.0, 0.0), (1.0, 5.0)])
        self.assertTrue(line0.intersects(line1))
        self.assertEqual(line0.intersections(line1), Multipoint([(0.5, 2.5)]))
        return

    def test_line_intersection_vertical(self):
        line0 = Line([(2.5, 2.5), (2.5, -2.5)])
        line1 = Line([(1.5, 2.5), (3.5, -2.5)])
        self.assertTrue(line0.intersects(line1))
        self.assertEqual(line0.intersections(line1), Multipoint([(2.5, 0.0)]))
        return

    def test_intersection_polygons(self):
        poly0 = Polygon([(0, 0), (2, 0), (3, 1), (2, 1), (2, 2), (1, 0)])
        poly1 = Polygon([(-1, -1), (1, -1), (1, 1), (-1, 1)])
        self.assertTrue(poly0.intersects(poly1))
        return

    def test_line_intersects_geographical1(self):
        line1 = Line([(-40.0, 36.0), (-38.0, 36.5)], crs=SphericalEarth)
        line2 = Line([(-39.0, 34.0), (-39.0, 37.5)], crs=SphericalEarth)
        self.assertTrue(line1.intersects(line2))
        return

    def test_line_intersects_geographical2(self):
        line1 = Line([(-40.0, 36.0), (-38.0, 36.5)], crs=SphericalEarth)
        line2 = Line([(-42.0, 34.0), (-41.0, 37.5)], crs=SphericalEarth)
        self.assertFalse(line1.intersects(line2))
        return

    def test_line_intersects_geographical3(self):
        # checks to make sure geodesics are handled
        line1 = Line([(-50.0, 70.0), (50.0, 70.0)], crs=SphericalEarth)
        line2 = Line([(0.0, 71.0), (1.0, 89.0)], crs=SphericalEarth)
        self.assertTrue(line1.intersects(line2))
        return

    def test_line_intersects_geographical4(self):
        # catches possible bugs in handling vertical segments on sweepline
        line1 = Line([(-50.0, 70.0), (50.0, 70.0)], crs=SphericalEarth)
        line2 = Line([(0.0, 71.0), (0.0, 89.0)], crs=SphericalEarth)
        self.assertTrue(line1.intersects(line2))
        return

    def test_poly_contains1(self):
        # trivial cases
        pt0 = Point((-0.5, 0.92))
        unitsquare = Polygon([(0.0,0.0), (1.0,0.0), (1.0,1.0), (0.0,1.0)])
        self.assertFalse(unitsquare.contains(pt0))

        pt1 = Point((0.125, 0.875))
        self.assertTrue(unitsquare.contains(pt1))

        x = np.arange(-4, 5)
        y = (x)**2
        line = Line([(x_,y_) for x_,y_ in zip(x, y)], crs=Cartesian)
        bbox = Polygon([(-2.5, 2.5), (2.5, 2.5), (2.5, -2.5), (-2.5, -2.5)],
                             crs=Cartesian)

        self.assertEqual(list(filter(bbox.contains, line)),
                         [Point((-1, 1)), Point((0, 0)), Point((1, 1))])
        return

    def test_poly_contains2(self):
        # test some hard cases
        diamond = Polygon([(0,0), (1,1), (2,0), (1, -1)])
        self.assertFalse(diamond.contains(Point((2, 1))))
        self.assertTrue(diamond.contains(Point((1, 0))))
        self.assertFalse(diamond.contains(Point((2.5, 0))))
        self.assertFalse(diamond.contains(Point((0, -1))))
        self.assertFalse(diamond.contains(Point((2, -1))))
        return

    def test_poly_contains3(self):
        # case where point is on an edge (should return true)
        square = Polygon([(0,0), (1,0), (1,1), (0,1)])
        self.assertTrue(square.contains(Point([0.5, 0])))
        self.assertTrue(square.contains(Point([0, 0.5])))
        return

    def test_poly_contains4(self):
        # hippie star
        theta = np.linspace(0, 2*np.pi, 361)[:-1]
        r = 10*np.sin(theta*8) + 15
        x = np.cos(theta) * r + 25
        y = np.sin(theta) * r + 25
        polygon = Polygon(zip(x, y))
        # causes naive cross-product methods to fail
        pt = Point((28.75, 25.625))
        self.assertTrue(polygon.contains(pt))
        return

    def test_poly_contains_polar(self):
        p = Polygon([(0, 80), (45, 80), (90, 80), (135, 80), (180, 80),
                     (225, 80), (270, 80), (315, 80)],
                    crs=SphericalEarth)
        self.assertTrue(p.contains(Point((45, 85), crs=SphericalEarth)))
        self.assertFalse(p.contains(Point((45, 75), crs=SphericalEarth)))
        return

    def test_within_distance(self):
        line = Line([(0,0), (1,1), (3,1)])
        pt = Point((1,1.5))
        self.assertTrue(line.within_distance(pt, 0.6))
        self.assertFalse(line.within_distance(pt, 0.4))
        return

    def test_multipoint_within_bbox(self):
        vertices = [(float(x),float(y)) for x in range(-10,11)
                                        for y in range(-10,11)]
        ans = [v for v in vertices if (-5.0<v[0]<5.0) and (-4.0<v[1]<6.0)]
        mp = Multipoint(vertices)
        sub = mp.within_bbox((-5.0, -4.0, 5.0, 6.0))
        self.assertEqual(sub, Multipoint(ans))
        return

    def test_multipoint_within_polygon(self):
        np.random.seed(42)
        x = (np.random.random(100) - 0.5) * 180.0
        y = (np.random.random(100) - 0.5) * 30.0
        xp = [-80, -50, 20, 35, 55, -45, -60]
        yp = [0, -10, -8, -17, 15, 18, 12]
        poly = Polygon(zip(xp, yp), crs=LonLatWGS84)
        mp = Multipoint(zip(x, y), crs=LonLatWGS84)

        subset = mp.within_polygon(poly)
        excluded = [pt for pt in mp if pt not in subset]
        self.assertTrue(all(poly.contains(pt) for pt in subset))
        self.assertFalse(any(poly.contains(pt) for pt in excluded))
        return

    def test_multiline_touching_line(self):
        np.random.seed(49)
        multiline = Multiline([10*np.random.rand(10, 2)
                               + np.random.randint(-50, 50, (1, 2)) for _ in range(50)])
        line = Line([(-30, -40), (11, -30), (10, 22), (-10, 50)])
        touching = multiline.touching(line)
        self.assertEqual(len(touching), 4)
        return

    def test_multipolygon_touching_line(self):
        np.random.seed(49)
        multipolygon = \
            Multipolygon([[np.array([[0,0],[10,0],[10,10],[0,10]])
                           + np.random.randint(-50, 50, (1, 2))]
                          for _ in range(50)])
        line = Line([(-40, -35), (-15, -30), (30, 5), (10, 32), (-15, 17)])
        touching = multipolygon.touching(line)
        self.assertEqual(len(touching), 10)
        return

    def test_multiline_touching_poly(self):
        np.random.seed(49)
        multiline = Multiline([10*np.random.rand(10, 2)
                               + np.random.randint(-50, 50, (1, 2)) for _ in range(50)])
        poly = Polygon([(-30, -40), (12, -30), (8, 22), (-10, 50)])
        touching = multiline.touching(poly)
        self.assertEqual(len(touching), 12)
        return

    def test_multipolygon_touching_poly(self):
        np.random.seed(49)
        multipolygon = \
            Multipolygon([[np.array([[0,0],[3,0],[3,3],[0,3]])
                           + np.random.randint(-50, 50, (1, 2))]
                          for _ in range(50)])
        poly = Polygon([(-30, -40), (12, -30), (8, 22), (-10, 50)])
        touching = multipolygon.touching(poly)
        self.assertEqual(len(touching), 14)
        return

    def test_multiline_within_poly(self):
        np.random.seed(49)
        multiline = Multiline([10*np.random.rand(10, 2)
                               + np.random.randint(-50, 50, (1, 2)) for _ in range(50)])
        poly = Polygon([(-30, -40), (12, -30), (8, 22), (-10, 50)])
        within = multiline.within(poly)
        self.assertEqual(len(within), 8)
        return

    def test_multipolygon_within_poly(self):
        np.random.seed(49)
        multipolygon = \
            Multipolygon([[np.array([[0,0],[3,0],[3,3],[0,3]])
                           + np.random.randint(-50, 50, (1, 2))]
                          for _ in range(50)])
        poly = Polygon([(-30, -40), (12, -30), (8, 22), (-10, 50)])
        within = multipolygon.within(poly)
        self.assertEqual(len(within), 8)
        return

if __name__ == "__main__":
    unittest.main()
