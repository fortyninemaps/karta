""" Unit tests for vector functions """

import unittest
import os
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import json
from test_helper import TESTDATA

import karta.vector as vector
import karta.vector._geojson as geojson
from karta.vector._geojson import GeoJSONReader, GeoJSONNamedCRS
from karta.vector.geometry import Point, Line, Polygon, Multipoint, Multiline, Multipolygon
from karta.crs import LonLatWGS84, WebMercator, Cartesian

class TestGeoJSONInput(unittest.TestCase):

    def test_point_read(self):
        with open(os.path.join(TESTDATA, 'geojson_input/point.json')) as f:
            reader = GeoJSONReader(f)
        res = reader.parse()
        self.assertEqual(res.coordinates, [100.0, 0.0])
        return

    def test_linestring_read(self):
        with open(os.path.join(TESTDATA, 'geojson_input/linestring.json')) as f:
            reader = GeoJSONReader(f)
        res = reader.parse()
        self.assertEqual(res.coordinates, [[100.0, 0.0], [101.0, 1.0]])
        return

    def test_polygon_read(self):
        with open(os.path.join(TESTDATA, 'geojson_input/polygon.json')) as f:
            reader = GeoJSONReader(f)
        res = reader.parse()
        self.assertEqual(res.coordinates,
            [[[100.0, 0.0], [101.0, 0.0], [101.0, 1.0], [100.0, 1.0], [100.0, 0.0]]])
        return

    def test_multipoint_read(self):
        with open(os.path.join(TESTDATA, 'geojson_input/multipoint.json')) as f:
            reader = GeoJSONReader(f)
        res = reader.parse()
        self.assertEqual(res.coordinates, [[100.0, 0.0], [101.0, 1.0]])
        return

    def test_multilinestring_read(self):
        with open(os.path.join(TESTDATA, 'geojson_input/multilinestring.json')) as f:
            reader = GeoJSONReader(f)
        res = reader.parse()
        self.assertEqual(res.coordinates, [[[100.0, 0.0], [101.0, 1.0]],
                                              [[102.0, 2.0], [103.0, 3.0]]])
        return

    def test_multipolygon_read(self):
        with open(os.path.join(TESTDATA, 'geojson_input/multipolygon.json')) as f:
            reader = GeoJSONReader(f)
        res = reader.parse()
        self.assertEqual(res.coordinates,
            [[[[102.0, 2.0], [103.0, 2.0], [103.0, 3.0], [102.0, 3.0], [102.0, 2.0]]],
             [[[100.0, 0.0], [101.0, 0.0], [101.0, 1.0], [100.0, 1.0], [100.0, 0.0]],
             [[100.2, 0.2], [100.8, 0.2], [100.8, 0.8], [100.2, 0.8], [100.2, 0.2]]]])
        return

    def test_featurecollection_read(self):
        path = os.path.join(TESTDATA, "geojson_input/featurecollection.json")
        with open(path) as f:
            reader = GeoJSONReader(f)
        fc = reader.parse()
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

    def setUp(self):
        self.default_crs = GeoJSONNamedCRS("urn:ogc:def:crs:OGC:1.3:CRS84")

    def test_geometrycollection2geometry(self):
        path = os.path.join(TESTDATA, "geojson_input/geometrycollection.json")
        geoms = vector.read_geojson(path)

        self.assertEqual(len(geoms), 2)
        self.assertTrue(isinstance(geoms[0], vector.Point))
        self.assertTrue(isinstance(geoms[1], vector.Line))
        return

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
                        properties={"prop0":"value0", "prop1":{"this":"that"}},
                        crs=LonLatWGS84)
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
        self.assertEqual(names, features[0].data.getfield("n"))
        return

    def test_read_with_crs(self):
        path = os.path.join(TESTDATA, "geojson_input/us-capitols.json")
        features = vector.read_geojson(path, crs=LonLatWGS84)
        for f in features:
            self.assertEqual(f.crs, LonLatWGS84)
        return

class TestGeoJSONOutput(unittest.TestCase):

    maxDiff = None

    def verifyJson(self, json1, json2):
        """ Verify that two JSON strings are equivalent """
        obj1 = json.loads(json1)
        obj2 = json.loads(json2)
        self.assertEqual(obj1, obj2)
        return

    def test_point_write_cartesian(self):
        p = Point((100.0, 0.0), crs=Cartesian)
        s = p.as_geojson(urn="urn:ogc:def:crs:EPSG::5806", force_wgs84=False)
        ans = """{"properties": {}, "geometry": {"coordinates": [100.0, 0.0], "crs": {"properties": {"name": "urn:ogc:def:crs:EPSG::5806"}, "type": "name"}, "type": "Point"}, "type": "Feature"}"""
        self.verifyJson(s, ans)
        return

    def test_point_write(self):
        p = Point((100.0, 0.0), crs=LonLatWGS84)
        s = p.as_geojson(urn="urn:ogc:def:crs:EPSG::5806")
        ans = """{"properties": {}, "geometry": {"coordinates": [100.0, 0.0], "crs": {"properties": {"name": "urn:ogc:def:crs:EPSG::5806"}, "type": "name"}, "type": "Point"}, "type": "Feature"}"""
        self.verifyJson(s, ans)
        return

    def test_line_write(self):
        p = Line([(100.0, 0.0), (101.0, 1.0)], crs=LonLatWGS84)
        s = p.as_geojson(urn="urn:ogc:def:crs:EPSG::5806")
        ans = """{ "type": "Feature", "properties": {}, "geometry": { "crs": { "type": "name", "properties": { "name": "urn:ogc:def:crs:EPSG::5806" } }, "coordinates": [ [ 100.0, 0.0 ], [ 101.0, 1.0 ] ], "type": "LineString" } }"""

        self.verifyJson(s, ans)
        return

    def test_polygon_write(self):
        p = Polygon([[100.0, 0.0], [101.0, 0.0], [101.0, 1.0],
                     [100.0, 1.0]], crs=LonLatWGS84)
        s = p.as_geojson(urn="urn:ogc:def:crs:EPSG::5806")
        ans = """{ "properties": {}, "geometry": { "type": "Polygon", "crs": { "type": "name", "properties": { "name": "urn:ogc:def:crs:EPSG::5806" } }, "coordinates": [ [ [ 100.0, 0.0 ], [ 101.0, 0.0 ], [ 101.0, 1.0 ], [ 100.0, 1.0 ], [ 100.0, 0.0 ] ] ] }, "type": "Feature" }"""

        self.verifyJson(s, ans)
        return

    def test_multiline_write(self):
        p = Multiline([[(100, 0), (101, 1)], [(102, 2), (103, 3)]], crs=LonLatWGS84)
        s = p.as_geojson(urn="urn:ogc:def:crs:EPSG::5806")
        ans = """{"type": "Feature", "properties": {}, "geometry" : { "crs": { "type": "name", "properties": { "name": "urn:ogc:def:crs:EPSG::5806" } }, "type": "MultiLineString", "coordinates": [ [ [100.0, 0.0], [101.0, 1.0] ], [ [102.0, 2.0], [103.0, 3.0] ] ] } }"""

        self.verifyJson(s, ans)

    def test_multipolygon_write(self):
        p = Multipolygon([[[(102, 2), (103, 2), (103, 3), (102, 3)]],
                          [[(100, 0), (101, 0), (101, 1), (100, 1)],
                           [(100.2, 0.2), (100.8, 0.2), (100.8, 0.8), (100.2, 0.8)]]],
                          crs=LonLatWGS84)
        s = p.as_geojson(urn="urn:ogc:def:crs:EPSG::5806")
        ans = """{"type": "Feature", "properties": {}, "geometry" : { "crs": { "type": "name", "properties": { "name": "urn:ogc:def:crs:EPSG::5806" } }, "type": "MultiPolygon",
    "coordinates": [
      [[[102.0, 2.0], [103.0, 2.0], [103.0, 3.0], [102.0, 3.0], [102.0, 2.0]]],
      [[[100.0, 0.0], [101.0, 0.0], [101.0, 1.0], [100.0, 1.0], [100.0, 0.0]],
       [[100.2, 0.2], [100.8, 0.2], [100.8, 0.8], [100.2, 0.8], [100.2, 0.2]]]
      ]
    } }"""

        self.verifyJson(s, ans)

    def test_write_reproject(self):
        # tests whether coordinates are correctly reprojected to WGS84 lon/lat
        p = Line([(1e6, 1e6), (1.2e6, 1.4e6)], crs=WebMercator)
        s = p.as_geojson()
        ans = """{ "type": "Feature", "properties": {},
            "geometry": {
                "crs": { "type": "name", "properties": { "name": "urn:ogc:def:crs:OGC:1.3:CRS84" } },
                "coordinates": [[8.983152841195214, 8.946573850543412],
                                [10.779783409434257, 12.476624651238847]],
                "type": "LineString" } }"""

        self.verifyJson(s, ans)
        return

    def test_write_string_data(self):
        capitols = Multipoint([(-112.1, 33.57), (-121.5, 38.57),
                        (-84.42, 33.76), (-86.15, 39.78), (-112.0, 46.6),
                        (-82.99, 39.98), (-77.48, 37.53), (-95.69, 39.04),
                        (-71.02, 42.33), (-96.68, 40.81)],
                        data = {"n": ["Phoenix, Arizona",
                                      "Sacramento, California",
                                      "Atlanta, Georgia",
                                      "Indianapolis, Indiana",
                                      "Helena, Montana", "Columbus, Ohio",
                                      "Richmond, Virginia", "Topeka, Kansas",
                                      "Boston, Massachusetts",
                                      "Lincoln, Nebraska"]},
                        crs=LonLatWGS84)
        s = capitols.as_geojson(urn="urn:ogc:def:crs:EPSG::5806")
        ans = """{ "properties": { "n": [ "Phoenix, Arizona", "Sacramento, California", "Atlanta, Georgia", "Indianapolis, Indiana", "Helena, Montana", "Columbus, Ohio", "Richmond, Virginia", "Topeka, Kansas", "Boston, Massachusetts", "Lincoln, Nebraska" ] }, "type": "Feature", "geometry": { "type": "MultiPoint", "crs": { "type": "name", "properties": { "name": "urn:ogc:def:crs:EPSG::5806" } }, "coordinates": [ [ -112.1, 33.57 ], [ -121.5, 38.57 ], [ -84.42, 33.76 ], [ -86.15, 39.78 ], [ -112.0, 46.6 ], [ -82.99, 39.98 ], [ -77.48, 37.53 ], [ -95.69, 39.04 ], [ -71.02, 42.33 ], [ -96.68, 40.81 ] ] } } """
        self.verifyJson(s, ans)
        return

    def test_write_data_crs(self):
        capitols = vector.geometry.multipart_from_singleparts(
                              [Point((-112.1, 33.57), crs=LonLatWGS84),
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
        self.assertTrue("crs" in s)
        self.assertTrue('"name": "urn:ogc:def:crs:OGC:1.3:CRS84"' in s)
        return

class GeoJSONSerializerTests(unittest.TestCase):

    def setUp(self):
        self.serializer = geojson.GeoJSONSerializer()
        return

    def test_serialize_point(self):
        pt = geojson.Point((44.0, 17.0), LonLatWGS84)
        s = self.serializer(pt)
        d = json.loads(s)
        self.assertEqual(tuple(pt.coordinates), tuple(d["coordinates"]))
        return

    def test_serialize_linestring(self):
        linestring = geojson.LineString([[44.0, 17.0], [43.0, 17.5], [-2.1, 4.0]],
                                        LonLatWGS84)
        s = self.serializer(linestring)
        d = json.loads(s)
        self.assertEqual(list(linestring.coordinates), list(d["coordinates"]))
        return

    def test_serialize_polygon(self):
        polygon = geojson.Polygon([[[44.0, 17.0], [43.0, 17.5], [-2.1, 4.0], [44.0, 17.0]],
                                   [[1.0, 1.0], [0.5, -0.5], [0.8, -0.7], [1.0, 1.0]]],
                                  LonLatWGS84)
        s = self.serializer(polygon)
        d = json.loads(s)

        self.assertEqual(list(polygon.coordinates), list(d["coordinates"]))
        return

    def test_serialize_multipoint(self):
        multipoint = geojson.MultiPoint([[44.0, 17.0], [43.0, 17.5], [-2.1, 4.0]],
                                        LonLatWGS84)
        s = self.serializer(multipoint)
        d = json.loads(s)
        self.assertEqual(list(multipoint.coordinates), list(d["coordinates"]))
        return

    def test_serialize_multilinestring(self):
        multilinestring = geojson.MultiLineString(
                            [[[44.0, 17.0], [43.0, 17.5], [-2.1, 4.0]],
                             [[49.0, -3.0], [48.0, -2.5], [2.9, -16.0]]],
                            LonLatWGS84)
        s = self.serializer(multilinestring)
        d = json.loads(s)
        self.assertEqual(list(multilinestring.coordinates), list(d["coordinates"]))
        return

    def test_serialize_multipolygon(self):
        multipolygon = geojson.MultiPolygon(
                            [[[44.0, 17.0], [43.0, 17.5], [-2.1, 4.0]],
                              [[1.0, 1.0], [0.5, -0.5], [0.8, [-0.7]]],
                             [[[49.0, -3.0], [48.0, -2.5], [2.9, -16.0]]]],
                            LonLatWGS84)
        s = self.serializer(multipolygon)
        d = json.loads(s)
        self.assertEqual(list(multipolygon.coordinates), list(d["coordinates"]))
        return

if __name__ == "__main__":
    unittest.main()
