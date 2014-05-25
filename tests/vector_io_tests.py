""" Unit tests for vector functions """

import unittest
import os
import numpy as np
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import json
from test_helper import md5sum, md5sum_file, TESTDATA, TESTDIR

import karta
import karta.vector as vector
import karta.crs as crs
import karta.vector.geojson as geojson
from karta.vector.geojson import GeoJSONReader
from karta.vector.guppy import Point, Multipoint, Line, Polygon
from karta.vector.metadata import Metadata


class TestGeoJSONInput(unittest.TestCase):

    def test_point_read(self):
        with open(os.path.join(TESTDATA, 'geojson_input/point.json')) as f:
            reader = GeoJSONReader(f)
        res = reader.pull_points()
        self.assertEqual(res[0].coordinates, [100.0, 0.0])
        return

    def test_linestring_read(self):
        with open(os.path.join(TESTDATA, 'geojson_input/linestring.json')) as f:
            reader = GeoJSONReader(f)
        res = reader.pull_lines()
        self.assertEqual(res[0].coordinates, [[100.0, 0.0], [101.0, 1.0]])
        return

    def test_polygon_read(self):
        with open(os.path.join(TESTDATA, 'geojson_input/polygon.json')) as f:
            reader = GeoJSONReader(f)
        res = reader.pull_polygons()
        self.assertEqual(res[0].coordinates,
            [[[100.0, 0.0], [101.0, 0.0], [101.0, 1.0], [100.0, 1.0], [100.0, 0.0]]])
        return

    def test_multipoint_read(self):
        with open(os.path.join(TESTDATA, 'geojson_input/multipoint.json')) as f:
            reader = GeoJSONReader(f)
        res = reader.pull_multipoints()
        self.assertEqual(res[0].coordinates, [[100.0, 0.0], [101.0, 1.0]])
        return

    def test_multilinestring_read(self):
        with open(os.path.join(TESTDATA, 'geojson_input/multilinestring.json')) as f:
            reader = GeoJSONReader(f)
        res = reader.pull_lines()
        self.assertEqual(res[0].coordinates, [[100.0, 0.0], [101.0, 1.0]])
        self.assertEqual(res[1].coordinates, [[102.0, 2.0], [103.0, 3.0]])
        return

    def test_multipolygon_read(self):
        with open(os.path.join(TESTDATA, 'geojson_input/multipolygon.json')) as f:
            reader = GeoJSONReader(f)
        res = reader.pull_polygons()
        self.assertEqual(res[0].coordinates, 
            [[[[102.0, 2.0], [103.0, 2.0], [103.0, 3.0], [102.0, 3.0], [102.0, 2.0]]],
             [[[100.0, 0.0], [101.0, 0.0], [101.0, 1.0], [100.0, 1.0], [100.0, 0.0]],
             [[100.2, 0.2], [100.8, 0.2], [100.8, 0.8], [100.2, 0.8], [100.2, 0.2]]]])
        return

    def test_featurecollection_read(self):
        path = os.path.join(TESTDATA, "geojson_input/featurecollection.json")
        with open(path) as f:
            reader = GeoJSONReader(f)
        features = reader.pull_features()
        self.assertTrue(isinstance(features[0][0][0], geojson.Point))
        self.assertEqual(features[0][0][0].coordinates, [102.0, 0.5])
        self.assertEqual(features[0][1], {"prop0": "value0"})

        self.assertTrue(isinstance(features[1][0][0], geojson.LineString))
        self.assertEqual(features[1][0][0].coordinates, [[102.0, 0.0],
                             [103.0, 1.0], [104.0, 0.0], [105.0, 1.0]])
        self.assertEqual(features[1][1], {"prop0": "value0", "prop1": 0.0})

        self.assertTrue(isinstance(features[2][0][0], geojson.Polygon))
        self.assertEqual(features[2][0][0].coordinates, [[[100.0, 0.0],
                [101.0, 0.0], [101.0, 1.0], [100.0, 1.0], [100.0, 0.0]]])
        self.assertEqual(features[2][1], {"prop0": "value0",
                                          "prop1": {"this": "that"}})
        return


#class TestGuppyGeoJSON(unittest.TestCase):
#    """ While the test cases TestGeoJSONInput and TestGeoJSONOutput test the
#    low level geojson module, this test case focuses on the bindings with guppy.
#    """
#
#    def test_featurecollection_read(self):
#        path = os.path.join(TESTDATA, "geojson_input/featurecollection.json")
#        features = vector.read_geojson_features(path)
#        print(features)


class TestGeoJSONOutput(unittest.TestCase):

    def verifyJson(self, json1, json2):
        """ Verify that two JSON strings are equivalent """
        obj1 = json.loads(json1)
        obj2 = json.loads(json2)
        self.assertEqual(obj1, obj2)
        return

    def asJsonBuffer(self, geo):
        s = StringIO()
        geo.to_geojson(s)
        s.seek(0)
        return s

    def test_point_write(self):
        p = Point((100.0, 0.0))
        s = self.asJsonBuffer(p)
        ans = """{ "crs": { "properties": { "name": "urn:ogc:def:crs:EPSG::5806" }, "type": "name" }, "coordinates": [ 100.0, 0.0 ], "type": "Point" }"""
        self.verifyJson(s.read(), ans)
        return

    def test_line_write(self):
        p = Line([(100.0, 0.0), (101.0, 1.0)])
        s = self.asJsonBuffer(p)
        ans = """{ "type": "Feature", "crs": { "type": "name", "properties": { "name": "urn:ogc:def:crs:EPSG::5806" } }, "properties": {}, "id": [ 0, 1 ], "geometry": { "coordinates": [ [ 100.0, 0.0 ], [ 101.0, 1.0 ] ], "type": "LineString" }, "bbox": [ [ 100.0, 101.0 ], [ 0.0, 1.0 ] ] }"""

        self.verifyJson(s.read(), ans)
        return

    def test_polygon_write(self):
        p = Polygon([[100.0, 0.0], [101.0, 0.0], [101.0, 1.0],
                            [100.0, 1.0], [100.0, 0.0]])
        s = self.asJsonBuffer(p)
        ans = """{ "bbox": [ [ 100.0, 101.0 ], [ 0.0, 1.0 ] ], "properties": {}, "id": [ 0, 1, 2, 3, 4 ], "crs": { "type": "name", "properties": { "name": "urn:ogc:def:crs:EPSG::5806" } }, "geometry": { "type": "Polygon", "coordinates": [ [ [ 100.0, 0.0 ], [ 101.0, 0.0 ], [ 101.0, 1.0 ], [ 100.0, 1.0 ], [ 100.0, 0.0 ] ] ] }, "type": "Feature" }"""
        self.verifyJson(s.read(), ans)
        return

    def test_write_string_data(self):
        capitols = Multipoint(
            [Point((-112.1, 33.57), properties={"n": "Phoenix, Arizona"}),
             Point((-121.5, 38.57), properties={"n": "Sacramento, California"}),
             Point((-84.42, 33.76), properties={"n": "Atlanta, Georgia"}),
             Point((-86.15, 39.78), properties={"n": "Indianapolis, Indiana"}),
             Point((-112.0, 46.6,) , properties={"n": "Helena, Montana"}),
             Point((-82.99, 39.98), properties={"n": "Columbus, Ohio"}),
             Point((-77.48, 37.53), properties={"n": "Richmond, Virginia"}),
             Point((-95.69, 39.04), properties={"n": "Topeka, Kansas"}),
             Point((-71.02, 42.33), properties={"n": "Boston, Massachusetts"}),
             Point((-96.68, 40.81), properties={"n": "Lincoln, Nebraska"})])
        s = capitols.as_geojson()
        ans = """{ "properties": { "n": [ "Phoenix, Arizona", "Sacramento, California", "Atlanta, Georgia", "Indianapolis, Indiana", "Helena, Montana", "Columbus, Ohio", "Richmond, Virginia", "Topeka, Kansas", "Boston, Massachusetts", "Lincoln, Nebraska" ] }, "bbox": [ [ -121.5, -71.02 ], [ 33.57, 46.6 ] ], "type": "Feature", "crs": { "type": "name", "properties": { "name": "urn:ogc:def:crs:EPSG::5806" } }, "geometry": { "type": "MultiPoint", "coordinates": [ [ -112.1, 33.57 ], [ -121.5, 38.57 ], [ -84.42, 33.76 ], [ -86.15, 39.78 ], [ -112.0, 46.6 ], [ -82.99, 39.98 ], [ -77.48, 37.53 ], [ -95.69, 39.04 ], [ -71.02, 42.33 ], [ -96.68, 40.81 ] ] }, "id": [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 ] } """
        self.verifyJson(s, ans)
        return

    def test_write_data_crs(self):
        capitols = Multipoint([Point((-112.1, 33.57), crs=crs.LONLAT),
                               Point((-121.5, 38.57), crs=crs.LONLAT),
                               Point((-84.42, 33.76), crs=crs.LONLAT),
                               Point((-86.15, 39.78), crs=crs.LONLAT),
                               Point((-112.0, 46.6), crs=crs.LONLAT),
                               Point((-82.99, 39.98), crs=crs.LONLAT),
                               Point((-77.48, 37.53), crs=crs.LONLAT),
                               Point((-95.69, 39.04), crs=crs.LONLAT),
                               Point((-71.02, 42.33), crs=crs.LONLAT),
                               Point((-96.68, 40.81), crs=crs.LONLAT),
                               Point((-97.51, 35.47), crs=crs.LONLAT),
                               Point((-134.2, 58.37), crs=crs.LONLAT),
                               Point((-100.3, 44.38), crs=crs.LONLAT)])
        s = capitols.as_geojson()
        #print(s)


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
        track = [Line([(np.random.random(), np.random.random())
                       for i in range(10)], properties={"name":"segment0"})]
        g = vector.gpx.GPX()
        g.add_track(track)
        expected = self.Track([self.Trkseg(
                        [self.Point(xy, {}, {}) for xy in track[0].vertices],
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

