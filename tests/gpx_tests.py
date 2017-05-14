import unittest
import os
from test_helper import TESTDATA

import numpy as np
import karta
import karta.vector as vector
from karta.vector.geometry import Point, Multipoint, Line, Polygon

class GPXTests(unittest.TestCase):

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
        track = Multipoint([(np.random.random(), np.random.random())
                      for i in range(10)], properties={"name":"segment0"})
        g = vector.gpx.GPX()
        g.add_track(track)
        expected = self.Track([self.Trkseg(
                        [self.Point(tuple(xy), {}, {}) for xy in track.vertices],
                        {"name":"segment0"}, {})], {}, {})
        self.assertEqual(g.tracks[0], expected)
        return

    def test_add_route(self):
        route = Multipoint([(np.random.random(), np.random.random())
                      for i in range(10)], properties={"name":"route0"})
        g = vector.gpx.GPX()
        g.add_route(route)
        expected = self.Route([self.Point(tuple(xy), {}, {}) for xy in route.vertices],
                              {"name":"route0"}, {})
        self.assertEqual(g.routes[0], expected)
        return

class ReadGPXTests(unittest.TestCase):

    def test_mtn_bike_trail(self):
        tracks = vector.read_gpx_tracks(os.path.join(TESTDATA, "gpx_input", "fishermans-trail.gpx"))
        track1 = tracks[0]
        seg1 = track1[0]
        self.assertEqual(seg1.bbox(), (-123.00702, 49.32947, -122.991408, 49.392751))

if __name__ == "__main__":
    unittest.main()
