""" Unit tests for vector functions """

import unittest
import numpy as np
import shapely.geometry

import karta.vector as vector
from karta.vector.geometry import Point, Multipoint, Line, Polygon

class TestGeoInterface(unittest.TestCase):

    def test_point_output(self):
        p = Point((4, 2))
        sp = shapely.geometry.shape(p.geomdict)
        self.assertEqual(sp.x, p.x)
        self.assertEqual(sp.y, p.y)
        return

    def test_multipoint_output(self):
        p = Multipoint([(4, 2), (3, 5), (3, 2), (7, 3)])
        sp = shapely.geometry.shape(p.geomdict)
        x, y = p.coordinates
        self.assertEqual(x, tuple([el.x for el in sp]))
        self.assertEqual(y, tuple([el.y for el in sp]))
        return

    def test_line_output(self):
        p = Line([(4, 2), (3, 5), (3, 2), (7, 3)])
        sp = shapely.geometry.shape(p.geomdict)
        x, y = p.coordinates
        sx, sy = sp.xy
        self.assertEqual(x, tuple(sx))
        self.assertEqual(y, tuple(sy))
        return

    def test_poly_output(self):
        p = Polygon([(4, 2), (3, 5), (3, 2), (7, 3)])
        sp = shapely.geometry.shape(p.geomdict)
        self.assertEqual(p.bbox, sp.bounds)
        return

    def test_point_input(self):
        sp = shapely.geometry.Point((3,4))
        p = vector.read.from_shape(sp)
        self.assertEqual(p.x, sp.x)
        self.assertEqual(p.y, sp.y)
        return

    def test_line_input(self):
        sp = shapely.geometry.LineString([(3,4), (6,2), (2,5)])
        p = vector.read.from_shape(sp)
        x, y = p.coordinates
        sx, sy = sp.xy
        self.assertEqual(x, tuple(sx))
        self.assertEqual(y, tuple(sy))
        return

    def test_poly_input(self):
        sp = shapely.geometry.Polygon([(4, 2), (3, 5), (3, 2), (7, 3)])
        p = vector.read.from_shape(sp)
        self.assertEqual(p.bbox, sp.bounds)
        return

    def test_multipoly_input(self):
        sp1 = shapely.geometry.Polygon([(4, 2), (3, 5), (3, 2), (7, 3)])
        sp2 = shapely.geometry.Polygon([(7, 3), (9, 7), (2, 7), (2, 0)])
        smp = shapely.geometry.MultiPolygon([sp1, sp2])
        p1, p2 = vector.read.from_shape(smp)
        self.assertEqual(p1.bbox, sp1.bounds)
        self.assertEqual(p2.bbox, sp2.bounds)
        return

    def test_multiline_input(self):
        sp1 = shapely.geometry.LineString([(4, 2), (3, 5), (3, 2), (7, 3)])
        sp2 = shapely.geometry.LineString([(7, 3), (9, 7), (2, 7), (2, 0)])
        smp = shapely.geometry.MultiLineString([sp1, sp2])
        p1, p2 = vector.read.from_shape(smp)
        x, y = p1.coordinates
        sx, sy = sp1.xy
        self.assertEqual(x, tuple(sx))
        self.assertEqual(y, tuple(sy))
        x, y = p2.coordinates
        sx, sy = sp2.xy
        self.assertEqual(x, tuple(sx))
        self.assertEqual(y, tuple(sy))
        return

    def test_feature_input(self):
        class Pointy(object):
            __geo_interface__ = {'type': 'Point', 'coordinates': (0.0, 0.0)}
        class Placemark(object):
            __geo_interface__ = {
                'type': 'Feature',
                'properties': {'name': 'Phoo'},
                'geometry': Pointy.__geo_interface__ }
        p = vector.read.from_shape(Placemark())
        self.assertEqual(p.properties["name"], "Phoo")
        self.assertEqual(p.vertex, (0.0, 0.0))
        return

class TestGPX(unittest.TestCase):

    def setUp(self):
        self.points = [vector.gpx.Point((np.random.random(), np.random.random()),
                                        {}, {}) for i in range(20)]
        self.segments = [vector.gpx.Trkseg(self.points, {}, {})]
        self.tracks = [vector.gpx.Track(self.segments, {}, {})]
        self.routes = [vector.gpx.Route(self.points, {}, {})]

        self.Point = vector.gpx.Point
        self.Trkseg = vector.gpx.Trkseg
        self.Track = vector.gpx.Track
        self.Route = vector.gpx.Route
        return

    def test_track_init(self):
        """ Test initiation of a GPX file containing a single track. """
        g = vector.gpx.GPX()
        for i, pt in enumerate(self.points):
            g.waypts.append(pt)
        g.tracks = self.tracks
        g.routes = self.routes
        return

    def test_add_waypoint(self):
        waypoint = Point((-80.0, 82.0),
                                properties = {"name": "ellesmere", "ele": 100})
        g = vector.gpx.GPX()
        g.add_waypoint(waypoint)
        expected = self.Point((-80.0, 82.0), {"name": "ellesmere",
                                              "ele": "100"}, {})
        self.assertEqual(g.waypts[0], expected)
        return

    def test_add_track(self):
        track = Line([(np.random.random(), np.random.random())
                      for i in range(10)], properties={"name":"segment0"})
        g = vector.gpx.GPX()
        g.add_track(track)
        expected = self.Track([self.Trkseg(
                        [self.Point(xy, {}, {}) for xy in track.vertices],
                        {"name":"segment0"}, {})], {}, {})
        self.assertEqual(g.tracks[0], expected)
        return

    def test_add_route(self):
        route = Line([(np.random.random(), np.random.random())
                      for i in range(10)], properties={"name":"route0"})
        g = vector.gpx.GPX()
        g.add_route(route)
        expected = self.Route([self.Point(xy, {}, {}) for xy in route.vertices],
                              {"name":"route0"}, {})
        self.assertEqual(g.routes[0], expected)
        return

if __name__ == "__main__":
    unittest.main()
