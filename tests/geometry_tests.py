""" Unit tests for vector functions """

from __future__ import division
import unittest
import math
import numpy as np

from karta.vector.geometry import Point, Multipoint, Line, Polygon
from karta.vector.geometry import affine_matrix
from karta.crs import Cartesian, SphericalEarth, LonLatWGS84, NSIDCNorth, Proj4CRS
from karta.crs import CRSError

class TestGeometry(unittest.TestCase):

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
        self.assertEqual(point.azimuth(other), 0.25*180)

        other = Point((0.0, 3.0))
        self.assertEqual(point.azimuth(other), 1.75*180)

        other = Point((0.0, 1.0))
        self.assertEqual(point.azimuth(other), 1.25*180)

        other = Point((2.0, 1.0))
        self.assertEqual(point.azimuth(other), 0.75*180)

        other = Point((1.0, 3.0))
        self.assertEqual(point.azimuth(other), 0.0)

        other = Point((1.0, 1.0))
        self.assertEqual(point.azimuth(other), 180.0)
        return

    def test_point_azimuth2(self):
        point = Point((5.0, 2.0))
        other = Point((5.0, 2.0))
        self.assertTrue(np.isnan(point.azimuth(other)))
        return

    def test_point_azimuth3(self):
        """ Verify with:

        printf "0 -1000000\n100000 -900000" | proj +proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +units=m +datum=WGS84 +no_defs -I -s | tr '\n' ' ' | invgeod +ellps=WGS84 -f "%.6f"
        """
        point = Point((0.0, -10e5), crs=NSIDCNorth)
        other = Point((1e5, -9e5), crs=NSIDCNorth)
        self.assertAlmostEqual(point.azimuth(other), 45.036973, places=6)
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

    def test_empty_multipoint(self):
        mp = Multipoint([], crs=LonLatWGS84)
        self.assertEqual(len(mp), 0)
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
        mp1[2] = Point((4.0, 5.0), data=["umiujaq"])
        self.assertEqual(mp1, mp2)
        return

    def test_multipoint_iterator(self):
        mp = Multipoint([(3.0, 3.0), (5.0, 1.0), (3.0, 1.0),
                         (4.0, 4.0), (0.0, 1.0)],
                         data=["rankin", "corbet", "arviat",
                               "severn", "churchill"])
        for i, pt in enumerate(mp):
            self.assertEqual(mp[i], pt)
        return

    def test_multipoint_get_data_fields(self):
        mp = Multipoint([(3.0, 3.0), (5.0, 1.0), (3.0, 1.0),
                         (4.0, 4.0), (0.0, 1.0)],
                         data={"location": ["rankin", "corbet", "arviat",
                                            "severn", "churchill"]})
        pt = mp[3]
        self.assertEqual(pt.data.fields, ("location",))
        self.assertEqual(pt.data[0], ("severn",))
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
        bbox = (1.0, 0.0, 9.0, 9.0)
        self.assertEqual(self.mp.bbox, bbox)
        return

    def test_multipoint_bbox_overlap(self):
        self.assertTrue(self.mp._bbox_overlap(self.poly))
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
        sub = mp.within_bbox((-5.0, -4.0, 5.0, 6.0))
        self.assertEqual(sub, Multipoint(ans))
        return

    def test_multipoint_convex_hull(self):
        vertices = [(953, 198), (986, 271), (937, 305), (934, 464), (967, 595),
                (965, 704), (800, 407), (782, 322), (863, 979), (637, 689),
                (254, 944), (330, 745), (363, 646), (27, 990), (127, 696),
                (286, 352), (436, 205), (88, 254), (187, 85)]
        mp = Multipoint(vertices)
        ch = mp.convex_hull()
        hull_vertices = [(187, 85), (953, 198), (986, 271), (965, 704), (863,
            979), (27, 990), (88, 254)]
        self.assertEqual(ch.vertices, hull_vertices)
        return

    def test_multipoint_convex_hull2(self):
        vertices = [(-158, 175), (-179, 230), (-404, -390), (259, -79), (32,
            144), (-59, 355), (402, 301), (239, 159), (-421, 172), (-482, 26),
            (2, -499), (134, -72), (-412, -12), (476, 235), (-412, 40), (-198,
                -256), (314, 331), (431, -492), (325, -415), (-400, -491)]
        mp = Multipoint(vertices)
        ch = mp.convex_hull()
        hull_vertices = [(2, -499), (431, -492), (476, 235), (402, 301), (314,
            331), (-59, 355), (-421, 172), (-482, 26), (-400, -491)]
        self.assertEqual(ch.vertices, hull_vertices)
        return

    def test_connected_multipoint_shortest_distance_to(self):
        line = Line([(0.0, 0.0), (2.0, 2.0), (5.0, 4.0)])
        dist = line.shortest_distance_to(Point((0.0, 2.0)))
        self.assertTrue(abs(dist - math.sqrt(2)) < 1e-10)
        return

    def test_connected_multipoint_shortest_distance_to2(self):
        line = Line([(127.0, -35.0), (132.0, -28.0), (142.0, -29.0)], crs=LonLatWGS84)
        dist = line.shortest_distance_to(Point((98.0, -7.0), crs=LonLatWGS84))
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
        line = Line([(-40, 0.0), (35, 0.0)], crs=LonLatWGS84)
        npt = line.nearest_on_boundary(Point((30.0, 80.0), crs=LonLatWGS84))
        self.assertPointAlmostEqual(npt, Point((30.0, 0.0), crs=LonLatWGS84))
        return

    # def test_connected_multipoint_nearest_on_boundary3(self):
    #     # This is the test that tends to break naive root finding schemes
    #     line = Line([(-40, 0.0), (35, 0.0)], crs=LonLatWGS84)
    #     npt = line.nearest_on_boundary(Point((30.0, 1e-8), crs=LonLatWGS84))
    #     self.assertPointAlmostEqual(npt, Point((30.0, 0.0), crs=LonLatWGS84))
    #     return

    # def test_connected_multipoint_nearest_on_boundary4(self):
    #     line = Line([(-20.0, 32.0), (-26.0, 43.0), (-38.0, 39.0)], crs=LonLatWGS84)
    #     npt = line.nearest_on_boundary(Point((-34.0, 52.0), crs=LonLatWGS84))
    #     self.assertPointAlmostEqual(npt, Point((-27.98347, 42.456316), crs=LonLatWGS84))
    #     return

    def test_line_append2d(self):
        ln0 = Line([(3.0, 3.0), (5.0, 1.0), (3.0, 1.0)])
        ln1 = Line([(3.0, 3.0), (5.0, 1.0), (3.0, 1.0),
                    (0.0, 1.0)])
        ln0.append(Point((0.0, 1.0)))
        self.assertEqual(ln0, ln1)
        return

    def test_line_append3d(self):
        ln0 = Line([(3.0, 3.0, 2.0), (5.0, 1.0, 0.0), (3.0, 1.0, 5.0)])
        ln1 = Line([(3.0, 3.0, 2.0), (5.0, 1.0, 0.0), (3.0, 1.0, 5.0),
                    (0.0, 1.0, 3.0)])
        ln0.append(Point((0.0, 1.0, 3.0)))
        self.assertEqual(ln0, ln1)
        return

    def test_line_extend(self):
        ln0a = Line([(3.0, 3.0, 2.0), (5.0, 1.0, 0.0), (3.0, 1.0, 5.0)])
        ln0b = Line([(4.0, 4.0, 6.0), (0.0, 1.0, 3.0)])
        ln1 = Line([(3.0, 3.0, 2.0), (5.0, 1.0, 0.0), (3.0, 1.0, 5.0),
                    (4.0, 4.0, 6.0), (0.0, 1.0, 3.0)])
        ln0a.extend(ln0b)
        self.assertEqual(ln0a, ln1)

    def test_pop(self):
        ln = Multipoint([(3.0, 3.0, 2.0), (5.0, 1.0, 0.0), (3.0, 1.0, 5.0),
                         (4.0, 4.0, 6.0), (0.0, 1.0, 3.0)],
                        data=["red", "green", "blue", "chartreuse", "aquamarine"])
        lnresult = Multipoint([(3.0, 3.0, 2.0), (5.0, 1.0, 0.0), (3.0, 1.0, 5.0),
                               (0.0, 1.0, 3.0)],
                              data=["red", "green", "blue", "aquamarine"])
        pt = ln.pop(3)
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
                     (120.0, 80.0), (150.0, 80.0), (180.0, 80.0), (-150.0, 80.0),
                     (-120.0, 80.0), (-90.0, 80.0), (-60.0, 80.0), (-30.0, 80.0)],
                    crs=SphericalEarth)
        self.assertTrue(p.ispolar())

        p = Polygon([(45.0, 30.0), (40.0, 25.0), (45.0, 20.0), (35.0, 25.0)],
                    crs=SphericalEarth)
        self.assertFalse(p.ispolar())


        p = Polygon([(45.0, 30.0), (40.0, 25.0), (45.0, 20.0), (35.0, 25.0)],
                    crs=Cartesian)
        self.assertRaises(CRSError, p.ispolar)
        return

    def test_poly_extents(self):
        self.assertEqual(self.poly.get_extents(), (0.0, 6.0, 1.0, 8.0))
        return

    def test_poly_extents_foreign_crs(self):
        x, y = zip(*self.poly.get_vertices(crs=NSIDCNorth))
        self.assertEqual(self.poly.get_extents(NSIDCNorth), (min(x), max(x), min(y), max(y)))
        return

    def test_poly_length(self):
        self.assertEqual(self.poly.length, 19.430647008220866)
        return

    def test_poly_contains1(self):
        # trivial case
        pt0 = Point((-0.5, 0.92))
        pt1 = Point((0.125, 0.875))
        self.assertFalse(self.unitsquare.contains(pt0))
        self.assertTrue(self.unitsquare.contains(pt1))
        return

    def test_poly_contains2(self):
        # trivial but more interesting case
        x = np.arange(-4, 5)
        y = (x)**2
        line = Line([(x_,y_) for x_,y_ in zip(x, y)], crs=Cartesian)
        bbox = Polygon([(-2.5, 2.5), (2.5, 2.5), (2.5, -2.5), (-2.5, -2.5)],
                             crs=Cartesian)

        self.assertEqual(list(filter(bbox.contains, line)),
                         [Point((-1, 1)), Point((0, 0)), Point((1, 1))])

    def test_poly_contains3(self):
        # test some hard cases
        diamond = Polygon([(0,0), (1,1), (2,0), (1, -1)])
        self.assertFalse(diamond.contains(Point((2, 1))))
        self.assertTrue(diamond.contains(Point((1, 0))))
        self.assertFalse(diamond.contains(Point((2.5, 0))))
        self.assertFalse(diamond.contains(Point((2, -1))))
        return

    def test_poly_contains4(self):
        # case where point is on an edge (should return true)
        square = Polygon([(0,0), (1,0), (1,1), (0,1)])
        pt = Point([0.5, 0])
        self.assertTrue(square.contains(pt))
        return

    def test_poly_contains5(self):
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

    def test_poly_centroid(self):
        poly = Polygon([(0,0), (1,0), (1,1), (0,1)], properties={"name": "features1"})
        c = poly.centroid
        self.assertEqual(c.x, 0.5)
        self.assertEqual(c.y, 0.5)
        self.assertEqual(c.properties, poly.properties)
        return

    def test_poly_centroid2(self):
        poly = Polygon([(0,0), (1,0), (2,0.5), (1,1), (0,1)], properties={"name": "features1"})
        c = poly.centroid
        self.assertAlmostEqual(c.x, 7/9)
        self.assertEqual(c.y, 0.5)
        self.assertEqual(c.properties, poly.properties)
        return

    def test_ringedpoly_perimeter(self):
        self.assertEqual(round(self.ringed_poly.perimeter, 3), 50.246)
        return

    def test_ringedpoly_area(self):
        self.assertEqual(self.ringed_poly.area, 100 - self.ring.area)
        return

    def test_area_compute_pi(self):
        r = np.linspace(0, 2*np.pi, 10000)
        x = np.cos(r)
        y = np.sin(r)
        kp = Polygon(zip(x,y))
        self.assertAlmostEqual(kp.area, np.pi, places=6)
        return

    def test_segments(self):
        v = self.vertices
        self.assertEqual([tuple(a.vertices) for a in self.line.segments],
                         [(v[i], v[i+1]) for i in range(len(self.vertices)-1)])
        return

    def test_within_distance(self):
        line = Line([(0,0), (1,1), (3,1)])
        pt = Point((1,1.5))
        self.assertTrue(line.within_distance(pt, 0.6))
        self.assertFalse(line.within_distance(pt, 0.4))
        return

    def test_walk_cartesian(self):
        start = Point((-3, -4), crs=Cartesian)
        dest = start.walk(5.0, math.atan(3.0/4.0), radians=True)
        self.assertAlmostEqual(dest.x, 0.0)
        self.assertAlmostEqual(dest.y, 0.0)
        return

    def test_walk(self):
        start = Point((-123.1, 49.25), crs=LonLatWGS84)
        dest = start.walk(1e5, 80.0)
        self.assertAlmostEqual(dest.x, -121.743196, places=6)
        self.assertAlmostEqual(dest.y, 49.398187, places=6)
        return

    def test_walk_albers(self):
        AlaskaAlbers = Proj4CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 "
                                "+x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs",
                                "+ellps=GRS80")
        start = Point((-2658638, 2443580), crs=AlaskaAlbers)
        dest = start.walk(4500, 195.0)
        self.assertAlmostEqual(dest.x, -2662670.889, places=3)
        self.assertAlmostEqual(dest.y, 2441551.155, places=3)
        return

    def test_subsection_cartesian(self):
        line = Line([(0.0, 0.0), (1.0, 2.0), (3.0, -2.0), (4.0, -1.0),
                     (4.0, 3.0), (3.0, 2.0)])
        points = line.subsection(20)
        ans = [Point(v) for v in [(0.0, 0.0),
                                  (0.318619234003536, 0.637238468007072),
                                  (0.637238468007072, 1.274476936014144),
                                  (0.9558577020106079, 1.9117154040212159),
                                  (1.274476936014144, 1.4510461279717122),
                                  (1.59309617001768, 0.8138076599646402),
                                  (1.911715404021216, 0.17656919195756826),
                                  (2.230334638024752, -0.4606692760495037),
                                  (2.5489538720282883, -1.0979077440565757),
                                  (2.867573106031824, -1.7351462120636478),
                                  (3.294395938694146, -1.7056040613058538),
                                  (3.7981771815888177, -1.2018228184111823),
                                  (4.0, -0.5729663008226373),
                                  (4.0, 0.13948796534818164),
                                  (4.0, 0.8519422315190006),
                                  (4.0, 1.5643964976898195),
                                  (4.0, 2.2768507638606383),
                                  (4.0, 2.989305030031457),
                                  (3.5037812428946715, 2.503781242894671),
                                  (3.0, 2.0)]]
        for a,b in zip(points, ans):
            self.assertPointAlmostEqual(a, b)
        return

    def test_subsection_lonlat(self):
        line = Line([(0, 40), (120, 40)], crs=LonLatWGS84)
        points = line.subsection(20)
        ans = [Point(v, crs=LonLatWGS84) for v in [(0, 40),
                                  (4.006549675732082, 43.200316625343305),
                                  (8.44359845345209, 46.2434129228378),
                                  (13.382442375999254, 49.09308515921458),
                                  (18.894149336762318, 51.705248417290484),
                                  (25.03918819127435, 54.027440893063556),
                                  (31.85052685770255, 55.99968253476488),
                                  (39.31083346558522, 57.55771841446013),
                                  (47.329401349484314, 58.6395037346357),
                                  (55.7308352362257, 59.194673757153645),
                                  (64.26916476377436, 59.19467375715364),
                                  (72.67059865051574, 58.639503734635674),
                                  (80.68916653441482, 57.557718414460105),
                                  (88.14947314229748, 55.999682534764844),
                                  (94.96081180872568, 54.02744089306352),
                                  (101.10585066323772, 51.705248417290456),
                                  (106.61755762400078, 49.09308515921457),
                                  (111.55640154654793, 46.24341292283779),
                                  (115.99345032426793, 43.2003166253433),
                                  (120, 40)]]
        for a,b in zip(points, ans):
            self.assertPointAlmostEqual(a, b)
        return

    def test_subsection_lonlat_precision(self):

        line = Line([(-20.247017, 79.683933), (-20.0993, 79.887917),
            (-19.13705, 80.048567), (-18.680467, 80.089333), (-17.451917,
                80.14405), (-16.913233, 80.02715), (-16.631367, 80.022933),
            (-16.194067, 80.0168), (-15.915983, 80.020267), (-15.7763,
                80.021283)], crs=LonLatWGS84)

        for n in range(2, 30):
            self.assertEqual(len(line.subsection(n)), n)
        return

class TestGeoInterface(unittest.TestCase):

    def test_point(self):
        pt = Point((1,2))
        self.assertEqual(pt.__geo_interface__, {"type":"Point", "coordinates":(1,2)})
        pt.shift((2,2))
        self.assertEqual(pt.__geo_interface__, {"type":"Point", "coordinates":(3,4)})

    def test_poly(self):
        x = np.arange(5)
        y = x**2
        poly = Polygon(list(zip(x, y)))
        self.assertEqual(poly.__geo_interface__,
                         {"type":"Polygon",
                          "bbox":(0, 0, 4, 16),
                          "coordinates": [list(zip(x, y))]})

    def test_line(self):
        x = np.arange(10)
        y = x**2
        line = Line(zip(x, y))
        self.assertEqual(line.__geo_interface__,
                         {"type":"LineString",
                          "bbox":(0, 0, 9, 81),
                          "coordinates": list(zip(x,y))})
        line = line[:5]
        self.assertEqual(line.__geo_interface__,
                         {"type":"LineString", 
                          "bbox":(0, 0, 4, 16),
                          "coordinates": list(zip(x[:5],y[:5]))})

class TestGeometryProj(unittest.TestCase):

    def setUp(self):
        self.vancouver = Point((-123.1, 49.25), crs=LonLatWGS84)
        self.ottawa = Point((-75.69, 45.42), crs=LonLatWGS84)
        self.whitehorse = Point((-135.05, 60.72), crs=LonLatWGS84)
        return

    def test_greatcircle(self):
        d1 = self.vancouver.distance(self.ottawa)
        d2 = self.vancouver.distance(self.whitehorse)
        d3 = self.whitehorse.distance(self.ottawa)
        d4 = self.whitehorse.distance(self.vancouver)
        self.assertTrue(abs(d1 - 3549030.70541) < 1e-5)
        self.assertTrue(abs(d2 - 1483327.53922) < 1e-5)
        self.assertTrue(abs(d3 - 4151366.88185) < 1e-5)
        self.assertTrue(abs(d4 - 1483327.53922) < 1e-5)
        return

    def test_azimuth_lonlat(self):
        """ Verify with
        echo "49.25dN 123.1dW 45.42dW 75.69dW" | invgeod +ellps=WGS84 -f "%.6f"
        """
        az1 = self.vancouver.azimuth(self.ottawa)
        az2 = self.vancouver.azimuth(self.whitehorse)
        self.assertAlmostEqual(az1, 78.483344, places=6)
        self.assertAlmostEqual(az2, -26.135827, places=6)
        return

    def test_walk_lonlat(self):
        start = Point((-132.14, 54.01), crs=LonLatWGS84)
        dest = start.walk(5440.0, 106.8)
        self.assertAlmostEqual(dest.x, -132.0605910876)
        self.assertAlmostEqual(dest.y, 53.99584742821)
        return


class TestGeometryOutput(unittest.TestCase):

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

class TestAffineTransforms(unittest.TestCase):

    def setUp(self):
        self.square = Polygon([(0, 0), (0, 1), (1, 1), (1, 0)])
        return

    def test_translate(self):
        M = affine_matrix(Multipoint([(0,0), (1,0), (0,1)]),
                          Multipoint([(1,0), (2,0), (1,1)]))
        Mans = np.array([[1, 0, 1], [0, 1, 0], [0, 0, 1]])
        self.assertTrue(np.allclose(M, Mans))

        translated_square = self.square.apply_affine_transform(M)
        ans = np.array([[1, 0], [1, 1], [2, 1], [2, 0]])
        self.assertTrue(np.allclose(translated_square.get_vertices(), ans))
        return

    def test_rotate(self):
        s2 = math.sqrt(0.5)
        M = affine_matrix(Multipoint([(0,0), (1,0), (0,1)]),
                          Multipoint([(0,0), (s2,s2), (-s2,s2)]))
        Mans = np.array([[s2, -s2, 0], [s2, s2, 0], [0, 0, 1]])
        self.assertTrue(np.allclose(M, Mans))

        translated_square = self.square.apply_affine_transform(M)
        ans = np.array([[0, 0], [-s2, s2], [0, 2*s2], [s2, s2]])
        self.assertTrue(np.allclose(translated_square.get_vertices(), ans))
        return

    def test_stretch(self):
        M = affine_matrix(Multipoint([(0,0), (1,0), (0,1)]),
                          Multipoint([(0,0), (2,0), (0,2)]))
        Mans = np.array([[2, 0, 0], [0, 2, 0], [0, 0, 1]])
        self.assertTrue(np.allclose(M, Mans))

        translated_square = self.square.apply_affine_transform(M)
        ans = np.array([[0, 0], [0, 2], [2, 2], [2, 0]])
        self.assertTrue(np.allclose(translated_square.get_vertices(), ans))
        return

class VectorCRSTests(unittest.TestCase):

    def test_vertices_in_crs(self):
        point = Point((-123.0, 49.0), crs=SphericalEarth)
        self.assertEqual(point.get_vertex(SphericalEarth),
                         (-123.0, 49.0))

    def test_vertices_in_crs2(self):
        point = Point((-123.0, 49.0), crs=SphericalEarth)
        BCAlbers = Proj4CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 "
                    "+lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 "
                    "+units=m +no_defs", "+ellps=GRS80")
        self.assertEqual(point.get_vertex(BCAlbers),
                         (1219731.770879303, 447290.49891930853))
        return

    def test_vertices_in_crs3(self):
        line = Line([(2.0, 34.0),
                     (2.15, 34.2),
                     (2.7, 34.1)], crs=SphericalEarth)
        UTM31N = Proj4CRS("+proj=utm +zone=31 +ellps=WGS84 "
                    "+datum=WGS84 +units=m +no_defs", "+ellps=WGS84")
        self.assertEqual(line.get_coordinate_lists(UTM31N),
                    ((407650.39665729366, 421687.71905896586, 472328.1095127584), 
                     (3762606.6598763773, 3784658.467084308, 3773284.485241791)))
        return

class MetadataAttributeTests(unittest.TestCase):

    def test_metadataattribute_str(self):
        g = Line(zip(range(5), range(5, 0, -1)),
                 data={"a": range(5), "b": [a**2 for a in range(-2, 3)]})
        self.assertEqual(g.d["a"], list(range(5)))
        self.assertEqual(g.d["b"], [a**2 for a in range(-2, 3)])
        return

    def test_metadataattribute_int(self):
        g = Line(zip(range(5), range(5, 0, -1)),
                 data={"a": range(5), "b": [a**2 for a in range(-2, 3)]})
        self.assertEqual(g.d[3], {"a": 3, "b": 1})
        return

    def test_None_data(self):
        g = Line(zip(range(5), range(5, 0, -1)))
        f = lambda a: getattr(g, a)
        self.assertRaises(KeyError, f, "d")
        return

if __name__ == "__main__":
    unittest.main()

