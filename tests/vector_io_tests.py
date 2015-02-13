""" Unit tests for vector functions """

import unittest
import os
import numpy as np
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import json
import shapely.geometry
from test_helper import TESTDATA

import karta.vector as vector
import karta.vector.geojson as geojson
from karta.vector.geojson import GeoJSONReader
from karta.vector.geometry import Point, Multipoint, Line, Polygon
from karta.crs import LonLatWGS84

class TestGeoInterface(unittest.TestCase):

    def test_point_conversion(self):
        p = Point((4, 2))
        shapely.geometry.shape(p)

    def test_multipoint_conversion(self):
        p = Multipoint([(4, 2), (3, 5), (3, 2), (7, 3)])
        shapely.geometry.shape(p)

    def test_line_conversion(self):
        p = Line([(4, 2), (3, 5), (3, 2), (7, 3)])
        shapely.geometry.shape(p)

    def test_poly_conversion(self):
        p = Polygon([(4, 2), (3, 5), (3, 2), (7, 3)])
        shapely.geometry.shape(p)

class TestGeoJSONInput(unittest.TestCase):

    def test_point_read(self):
        with open(os.path.join(TESTDATA, 'geojson_input/point.json')) as f:
            reader = GeoJSONReader(f)
        res = reader.items()
        self.assertEqual(res[0].coordinates, [100.0, 0.0])
        return

    def test_linestring_read(self):
        with open(os.path.join(TESTDATA, 'geojson_input/linestring.json')) as f:
            reader = GeoJSONReader(f)
        res = reader.items()
        self.assertEqual(res[0].coordinates, [[100.0, 0.0], [101.0, 1.0]])
        return

    def test_polygon_read(self):
        with open(os.path.join(TESTDATA, 'geojson_input/polygon.json')) as f:
            reader = GeoJSONReader(f)
        res = reader.items()
        self.assertEqual(res[0].coordinates,
            [[[100.0, 0.0], [101.0, 0.0], [101.0, 1.0], [100.0, 1.0], [100.0, 0.0]]])
        return

    def test_multipoint_read(self):
        with open(os.path.join(TESTDATA, 'geojson_input/multipoint.json')) as f:
            reader = GeoJSONReader(f)
        res = reader.items()
        self.assertEqual(res[0].coordinates, [[100.0, 0.0], [101.0, 1.0]])
        return

    def test_multilinestring_read(self):
        with open(os.path.join(TESTDATA, 'geojson_input/multilinestring.json')) as f:
            reader = GeoJSONReader(f)
        res = reader.items()
        self.assertEqual(res[0].coordinates, [[[100.0, 0.0], [101.0, 1.0]],
                                              [[102.0, 2.0], [103.0, 3.0]]])
        return

    def test_multipolygon_read(self):
        with open(os.path.join(TESTDATA, 'geojson_input/multipolygon.json')) as f:
            reader = GeoJSONReader(f)
        res = reader.items()
        self.assertEqual(res[0].coordinates, 
            [[[[102.0, 2.0], [103.0, 2.0], [103.0, 3.0], [102.0, 3.0], [102.0, 2.0]]],
             [[[100.0, 0.0], [101.0, 0.0], [101.0, 1.0], [100.0, 1.0], [100.0, 0.0]],
             [[100.2, 0.2], [100.8, 0.2], [100.8, 0.8], [100.2, 0.8], [100.2, 0.2]]]])
        return

    def test_featurecollection_read(self):
        path = os.path.join(TESTDATA, "geojson_input/featurecollection.json")
        with open(path) as f:
            reader = GeoJSONReader(f)
        fc = reader.items()[0]
        self.assertTrue(isinstance(fc.features[0].geometry, geojson.Point))
        self.assertEqual(fc.features[0].geometry.coordinates, [102.0, 0.5])
        self.assertEqual(fc.features[0].properties["scalar"], {"prop0": "value0"})

        self.assertTrue(isinstance(fc.features[1].geometry, geojson.LineString))
        self.assertEqual(fc.features[1].geometry.coordinates,
                        [[102.0, 0.0], [103.0, 1.0], [104.0, 0.0], [105.0, 1.0]])
        self.assertEqual(fc.features[1].properties["scalar"],
                        {"prop0": "value0", "prop1": 0.0})

        self.assertTrue(isinstance(fc.features[2].geometry, geojson.Polygon))
        self.assertEqual(fc.features[2].geometry.coordinates,
                        [[[100.0, 0.0], [101.0, 0.0], [101.0, 1.0],
                          [100.0, 1.0], [100.0, 0.0]]])
        self.assertEqual(fc.features[2].properties["scalar"],
                        {"prop0": "value0", "prop1": {"this": "that"}})
        return


class TestGeoJSON(unittest.TestCase):
    """ While the test cases TestGeoJSONInput and TestGeoJSONOutput test the
    low level geojson module, this test case focuses on the bindings with geometry.
    """

    def test_featurecollection2geometry(self):
        path = os.path.join(TESTDATA, "geojson_input/featurecollection.json")
        features = vector.read_geojson(path)

        ans0 = Point((102.0, 0.5), properties={"prop0":"value0"}, crs=LonLatWGS84)
        self.assertEqual(features[0], ans0)

        ans1 = Line([(102.0, 0.0), (103.0, 1.0), (104.0, 0.0), (105.0, 1.0)],
                    properties={"prop0":"value0", "prop1":0.0}, crs=LonLatWGS84)
        self.assertEqual(features[1], ans1)

        ans2 = Polygon([(100.0, 0.0), (101.0, 0.0), (101.0, 1.0), (100.0, 1.0),
                        (100.0, 0.0)],
                        properties={"prop0":"value0", "prop1":{"this":"that"}})
        self.assertEqual(features[2], ans2)
        return

    def test_read_capitols(self):
        path = os.path.join(TESTDATA, "geojson_input/us-capitols.json")
        features = vector.read_geojson(path)
        names = ['Phoenix, Arizona, United States', 'Sacramento, California, United States', 
                 'Atlanta, Georgia, United States', 'Indianapolis, Indiana, United States', 
                 'Helena, Montana, United States', 'Columbus, Ohio, United States', 
                 'Richmond, Virginia, United States', 'Topeka, Kansas, United States', 
                 'Boston, Massachusetts, United States', 'Lincoln, Nebraska, United States', 
                 'Oklahoma City, Oklahoma, United States', 'Juneau, Alaska, United States', 
                 'Pierre, South Dakota, United States', 'Honolulu, Hawaii, United States', 
                 'Montgomery, Alabama, United States',
                 'Little Rock, Arkansas, United States', 'Denver, Colorado, United States', 
                 'Hartford, Connecticut, United States', 'Dover, Delaware, United States', 
                 'Washington, District of Columbia, United States', 
                 'Tallahassee, Florida, United States', 'Boise, Idaho, United States', 
                 'Springfield, Illinois, United States', 'Des Moines, Iowa, United States', 
                 'Frankfort, Kentucky, United States', 
                 'Baton Rouge, Louisiana, United States', 'Augusta, Maine, United States', 
                 'Annapolis, Maryland, United States', 'Lansing, Michigan, United States', 
                 'Saint Paul, Minnesota, United States', 
                 'Jackson, Mississippi, United States', 
                 'Jefferson City, Missouri, United States', 
                 'Carson City, Nevada, United States', 
                 'Concord, New Hampshire, United States', 
                 'Trenton, New Jersey, United States', 
                 'Santa Fe, New Mexico, United States', 'Albany, New York, United States', 
                 'Raleigh, North Carolina, United States', 
                 'Bismarck, North Dakota, United States', 'Salem, Oregon, United States', 
                 'Harrisburg, Pennsylvania, United States', 
                 'Providence, Rhode Island, United States', 
                 'Columbia, South Carolina, United States', 
                 'Nashville, Tennessee, United States', 
                 'Austin, Texas, United States', 'Salt Lake City, Utah, United States', 
                 'Montpelier, Vermont, United States', 'Olympia, Washington, United States', 
                 'Charleston, West Virginia, United States', 
                 'Madison, Wisconsin, United States', 'Cheyenne, Wyoming, United States']
        self.assertEqual(names, features[0].data["n"])
        return


class TestGeoJSONOutput(unittest.TestCase):

    def verifyJson(self, json1, json2):
        """ Verify that two JSON strings are equivalent """
        obj1 = json.loads(json1)
        obj2 = json.loads(json2)
        self.assertEqual(obj1, obj2)
        return

    def asJsonBuffer(self, geo, **kw):
        s = StringIO()
        geo.to_geojson(s, **kw)
        s.seek(0)
        return s

    def test_point_write(self):
        p = Point((100.0, 0.0))
        s = self.asJsonBuffer(p, urn="urn:ogc:def:crs:EPSG::5806")
        ans = """{ "crs": { "properties": { "name": "urn:ogc:def:crs:EPSG::5806" }, "type": "name" }, "coordinates": [ 100.0, 0.0 ], "type": "Point" }"""
        self.verifyJson(s.read(), ans)
        return

    def test_line_write(self):
        p = Line([(100.0, 0.0), (101.0, 1.0)])
        s = self.asJsonBuffer(p, urn="urn:ogc:def:crs:EPSG::5806")
        ans = """{ "type": "Feature", "crs": { "type": "name", "properties": { "name": "urn:ogc:def:crs:EPSG::5806" } }, "properties": {}, "id": [ 0, 1 ], "geometry": { "coordinates": [ [ 100.0, 0.0 ], [ 101.0, 1.0 ] ], "type": "LineString" }, "bbox": [ [ 100.0, 101.0 ], [ 0.0, 1.0 ] ] }"""

        self.verifyJson(s.read(), ans)
        return

    def test_polygon_write(self):
        p = Polygon([[100.0, 0.0], [101.0, 0.0], [101.0, 1.0],
                            [100.0, 1.0], [100.0, 0.0]])
        s = self.asJsonBuffer(p, urn="urn:ogc:def:crs:EPSG::5806")
        ans = """{ "bbox": [ [ 100.0, 101.0 ], [ 0.0, 1.0 ] ], "properties": {}, "id": [ 0, 1, 2, 3, 4 ], "crs": { "type": "name", "properties": { "name": "urn:ogc:def:crs:EPSG::5806" } }, "geometry": { "type": "Polygon", "coordinates": [ [ [ 100.0, 0.0 ], [ 101.0, 0.0 ], [ 101.0, 1.0 ], [ 100.0, 1.0 ], [ 100.0, 0.0 ] ] ] }, "type": "Feature" }"""
        self.verifyJson(s.read(), ans)
        return

    def test_write_string_data(self):
        capitols = Multipoint(
            [Point((-112.1, 33.57), data={"n": "Phoenix, Arizona"}),
             Point((-121.5, 38.57), data={"n": "Sacramento, California"}),
             Point((-84.42, 33.76), data={"n": "Atlanta, Georgia"}),
             Point((-86.15, 39.78), data={"n": "Indianapolis, Indiana"}),
             Point((-112.0, 46.6,) , data={"n": "Helena, Montana"}),
             Point((-82.99, 39.98), data={"n": "Columbus, Ohio"}),
             Point((-77.48, 37.53), data={"n": "Richmond, Virginia"}),
             Point((-95.69, 39.04), data={"n": "Topeka, Kansas"}),
             Point((-71.02, 42.33), data={"n": "Boston, Massachusetts"}),
             Point((-96.68, 40.81), data={"n": "Lincoln, Nebraska"})])
        s = capitols.as_geojson(urn="urn:ogc:def:crs:EPSG::5806")
        ans = """{ "properties": { "n": [ "Phoenix, Arizona", "Sacramento, California", "Atlanta, Georgia", "Indianapolis, Indiana", "Helena, Montana", "Columbus, Ohio", "Richmond, Virginia", "Topeka, Kansas", "Boston, Massachusetts", "Lincoln, Nebraska" ] }, "bbox": [ [ -121.5, -71.02 ], [ 33.57, 46.6 ] ], "type": "Feature", "crs": { "type": "name", "properties": { "name": "urn:ogc:def:crs:EPSG::5806" } }, "geometry": { "type": "MultiPoint", "coordinates": [ [ -112.1, 33.57 ], [ -121.5, 38.57 ], [ -84.42, 33.76 ], [ -86.15, 39.78 ], [ -112.0, 46.6 ], [ -82.99, 39.98 ], [ -77.48, 37.53 ], [ -95.69, 39.04 ], [ -71.02, 42.33 ], [ -96.68, 40.81 ] ] }, "id": [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 ] } """
        self.verifyJson(s, ans)
        return

    def test_write_data_crs(self):
        capitols = Multipoint([Point((-112.1, 33.57), crs=LonLatWGS84),
                               Point((-121.5, 38.57), crs=LonLatWGS84),
                               Point((-84.42, 33.76), crs=LonLatWGS84),
                               Point((-86.15, 39.78), crs=LonLatWGS84),
                               Point((-112.0, 46.6), crs=LonLatWGS84),
                               Point((-82.99, 39.98), crs=LonLatWGS84),
                               Point((-77.48, 37.53), crs=LonLatWGS84),
                               Point((-95.69, 39.04), crs=LonLatWGS84),
                               Point((-71.02, 42.33), crs=LonLatWGS84),
                               Point((-96.68, 40.81), crs=LonLatWGS84),
                               Point((-97.51, 35.47), crs=LonLatWGS84),
                               Point((-134.2, 58.37), crs=LonLatWGS84),
                               Point((-100.3, 44.38), crs=LonLatWGS84)])
        s = capitols.as_geojson()


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

