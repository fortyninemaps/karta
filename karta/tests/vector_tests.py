""" Unit tests for vector functions """

import unittest
import os
import math
import numpy as np
import karta.vector as vector
from karta.vector.geojson import GeoJSONReader
from test_helper import md5sum, TESTDATA, CURDIR

class TestGuppy(unittest.TestCase):

    def setUp(self):
        self.point = vector.guppy.Point((1.0, 2.0, 3.0),
                                        data={"color":(43,67,10)},
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

        self.mp = vector.guppy.Multipoint(self.vertices, data=self.data)
        self.line = vector.guppy.Line(self.vertices, data=self.data)
        self.poly = vector.guppy.Polygon([(0.0, 8.0), (0.0, 5.0), (6.0, 1.0)])
        self.ring = vector.guppy.Polygon([(2.0, 2.0), (4.0, 2.0), (3.0, 6.0)])
        self.ringed_poly = vector.guppy.Polygon([(0.0, 0.0), (10, 0.0),
                                                 (10.0, 10.0), (0.0, 10.0)],
                                                subs=[self.ring])
        self.unitsquare = vector.guppy.Polygon([(0.0,0.0), (1.0,0.0), (1.0,1.0),
                                                (0.0,1.0)])
        return


    def test_point_coordsxy(self):
        self.assertEqual(self.point.coordsxy(), (1.0, 2.0))
        self.assertEqual(self.point[0], 1.0)
        self.assertEqual(self.point[1], 2.0)
        return

    def test_point_azimuth(self):
        other = vector.guppy.Point((7.0, 8.0))
        self.assertEqual(self.point.azimuth(other), 225.0 / 180.0 * np.pi)
        return

    def test_point_shift(self):
        point = vector.guppy.Point((-3.0, 5.0, 2.5))
        point.shift((4.0, -3.0, 0.5))
        self.assertEqual(self.point, point)
        return

    def test_nearest_to(self):
        self.assertEqual(self.mp.nearest_to(self.point), self.mp[12])
        return

    def test_multipoint_getset(self):
        self.assertEqual(self.mp[0], vector.guppy.Point(self.vertices[0],
                                                        data=self.mp.data[0],
                                                        properties=self.mp.properties))
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
        vertices = [(2.0, 9.0, 9.0),
                    (4.0, 1.0, 9.0),
                    (4.0, 1.0, 5.0),
                    (2.0, 8.0, 0.0),
                    (9.0, 8.0, 4.0),
                    (1.0, 4.0, 6.0),
                    (7.0, 3.0, 4.0),
                    (2.0, 5.0, 3.0),
                    (1.0, 6.0, 6.0),
                    (8.0, 1.0, 0.0),
                    (5.0, 5.0, 1.0),
                    (4.0, 5.0, 7.0),
                    (3.0, 3.0, 5.0),
                    (9.0, 0.0, 9.0),
                    (6.0, 3.0, 8.0),
                    (4.0, 5.0, 7.0),
                    (9.0, 9.0, 4.0),
                    (1.0, 4.0, 7.0),
                    (1.0, 7.0, 8.0),
                    (9.0, 1.0, 6.0)]

        data0 = [99.0, 2.0, 60.0, 75.0, 71.0, 34.0, 1.0, 49.0, 4.0, 36.0, 47.0,
                 58.0, 65.0, 72.0, 4.0, 27.0, 52.0, 37.0, 95.0, 17.0]

        data1 = [54.0, 40.0, 77.0, 18.0, 84.0, 91.0, 61.0, 92.0, 19.0, 42.0,
                 50.0, 25.0, 11.0, 80.0, 59.0, 56.0, 32.0, 8.0, 88.0, 76.0]

        L2 = vector.guppy.Multipoint(vertices, data={'d0':data0, 'd1':data1})
        return

    def test_connected_multipoint_distance_to(self):
        line = vector.guppy.Line([(0.0, 0.0), (2.0, 2.0), (5.0, 4.0)])
        dist = line.distance_to(vector.guppy.Point((0.0, 2.0)))
        self.assertTrue(abs(dist - math.sqrt(2)) < 1e-10)
        return

    def test_connected_multipoint_nearest_on_boundary(self):
        line = vector.guppy.Line([(0.0, 0.0), (2.0, 2.0), (5.0, 4.0)])
        npt = line.nearest_on_boundary(vector.guppy.Point((0.0, 2.0)))
        self.assertEqual(npt, vector.guppy.Point((1.0, 1.0)))
        return

    def test_line_intersection(self):
        line0 = vector.guppy.Line([(0.0, 0.0), (3.0, 3.0)])
        line1 = vector.guppy.Line([(0.0, 3.0), (3.0, 0.0)])
        self.assertTrue(line0.intersects(line1))
        self.assertEqual(line0.intersections(line1), [(1.5, 1.5)])
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
        pt0 = vector.guppy.Point((-0.5, 0.92))
        pt1 = vector.guppy.Point((0.125, 0.875))
        self.assertFalse(self.unitsquare.contains(pt0))
        self.assertTrue(self.unitsquare.contains(pt1))
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


class TestGuppyProj(unittest.TestCase):

    def setUp(self):
        self.vancouver = vector.guppy.Point((-123.1, 49.25), crs=vector.guppy.LONLAT)
        self.ottawa = vector.guppy.Point((-75.69, 45.42), crs=vector.guppy.LONLAT)
        self.whitehorse = vector.guppy.Point((-135.05, 60.72), crs=vector.guppy.LONLAT)
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

class TestGuppyOutput(unittest.TestCase):

    def setUp(self):

        if not os.path.isdir('data'):
            os.mkdir('data')

        vertices = [(2.0, 9.0, 9.0),
                    (4.0, 1.0, 9.0),
                    (4.0, 1.0, 5.0),
                    (2.0, 8.0, 0.0),
                    (9.0, 8.0, 4.0),
                    (1.0, 4.0, 6.0),
                    (7.0, 3.0, 4.0),
                    (2.0, 5.0, 3.0),
                    (1.0, 6.0, 6.0),
                    (8.0, 1.0, 0.0),
                    (5.0, 5.0, 1.0),
                    (4.0, 5.0, 7.0),
                    (3.0, 3.0, 5.0),
                    (9.0, 0.0, 9.0),
                    (6.0, 3.0, 8.0),
                    (4.0, 5.0, 7.0),
                    (9.0, 9.0, 4.0),
                    (1.0, 4.0, 7.0),
                    (1.0, 7.0, 8.0),
                    (9.0, 1.0, 6.0)]

        data0 = [99.0, 2.0, 60.0, 75.0, 71.0, 34.0, 1.0, 49.0, 4.0, 36.0, 47.0,
                 58.0, 65.0, 72.0, 4.0, 27.0, 52.0, 37.0, 95.0, 17.0]

        data1 = [54.0, 40.0, 77.0, 18.0, 84.0, 91.0, 61.0, 92.0, 19.0, 42.0,
                 50.0, 25.0, 11.0, 80.0, 59.0, 56.0, 32.0, 8.0, 88.0, 76.0]
        self.mp = vector.guppy.Multipoint(vertices,
                                          data={'d0':data0, 'd1':data1})

    def test_mp2vtp(self):
        # Test VTK output for a Multipoint
        with open(os.path.join(CURDIR, 'data', 'testmp2vtp.vtp'), 'w') as f:
            self.mp.to_vtk(f)
        self.assertEqual(md5sum('data/testmp2vtp.vtp'),
                         md5sum(os.path.join(TESTDATA, 'testmp2vtp.vtp')))
        return

    def test_geojson(self):
        # Test GeoJSON output for a Multipoint
        with open('data/testgeojson.json', 'w') as f:
            self.mp.to_geojson(f)
        self.assertEqual(md5sum('data/testgeojson.json'),
                         md5sum(os.path.join(TESTDATA, 'testgeojson.json')))
        return

class TestMetadata(unittest.TestCase):

    def setUp(self):
        self.onefield = vector.metadata.Metadata(data=np.arange(200), singleton=False)
        self.multifield = vector.metadata.Metadata(
                data={"a":np.arange(200), "b":np.arange(200,400), "c":np.arange(200)**2}, singleton=False)
        return

    def test_indexing(self):
        self.assertEqual(self.onefield[10], 10)
        self.assertEqual(self.multifield[10], {"a":10, "b": 210, "c": 100})
        return

    def test_slicing(self):
        self.assertTrue(np.all(self.onefield[5:15] == np.arange(5, 15)))
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

class TestGeoJSONInput(unittest.TestCase):

    def test_point_read(self):
        with open(os.path.join(TESTDATA, 'point.geojson')) as f:
            reader = GeoJSONReader(f)
        res = reader.pull_points()
        self.assertEqual(res[0].coordinates, [100.0, 0.0])
        return

    def test_linestring_read(self):
        with open(os.path.join(TESTDATA, 'linestring.geojson')) as f:
            reader = GeoJSONReader(f)
        res = reader.pull_lines()
        self.assertEqual(res[0].coordinates, [[100.0, 0.0], [101.0, 1.0]])
        return

    def test_polygon_read(self):
        with open(os.path.join(TESTDATA, 'polygon.geojson')) as f:
            reader = GeoJSONReader(f)
        res = reader.pull_polygons()
        self.assertEqual(res[0].coordinates, [[[100.0, 0.0], [101.0, 0.0],
                                               [101.0, 1.0], [100.0, 1.0],
                                               [100.0, 0.0]]])
        return

    def test_multipoint_read(self):
        with open(os.path.join(TESTDATA, 'multipoint.geojson')) as f:
            reader = GeoJSONReader(f)
        res = reader.pull_multipoints()
        self.assertEqual(res[0].coordinates, [[100.0, 0.0], [101.0, 1.0]])
        return

    def test_multilinestring_read(self):
        with open(os.path.join(TESTDATA, 'multilinestring.geojson')) as f:
            reader = GeoJSONReader(f)
        res = reader.pull_lines()
        self.assertEqual(res[0].coordinates, [[100.0, 0.0], [101.0, 1.0]])
        self.assertEqual(res[1].coordinates, [[102.0, 2.0], [103.0, 3.0]])
        return

    # This test will fail until holes are implemented
    def test_multipolygon_read(self):
        with open(os.path.join(TESTDATA, 'multipolygon.geojson')) as f:
            reader = GeoJSONReader(f)
        res = reader.pull_polygons()
        self.assertEqual(res[0].coordinates, [[[[102.0, 2.0], [103.0, 2.0],
                                                [103.0, 3.0], [102.0, 3.0],
                                                [102.0, 2.0]]],
                                              [[[100.0, 0.0], [101.0, 0.0],
                                                [101.0, 1.0], [100.0, 1.0],
                                                [100.0, 0.0]],
                                               [[100.2, 0.2], [100.8, 0.2],
                                                [100.8, 0.8], [100.2, 0.8],
                                                [100.2, 0.2]]]])
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
        """ Test initiation and writing of a GPX file containing a single
        track. """
        g = vector.gpx.GPX()
        for i, pt in enumerate(self.points):
            g.waypts.append(pt)
        g.tracks = self.tracks
        g.routes = self.routes
        g.writefile("test.gpx")
        return

    def test_add_waypoint(self):
        waypoint = vector.guppy.Point((-80.0, 82.0), properties = {"name":
                                                                   "ellesmere",
                                                                   "ele":
                                                                   100})
        g = vector.gpx.GPX()
        g.add_waypoint(waypoint)
        expected = self.Point((-80.0, 82.0), {"name": "ellesmere",
                                              "ele": "100"}, {})
        self.assertEqual(g.waypts[0], expected)
        return

    def test_add_track(self):
        track = [vector.guppy.Line([(np.random.random(), np.random.random())
                           for i in range(10)], properties={"name":"segment0"})]
        g = vector.gpx.GPX()
        g.add_track(track)
        expected = self.Track([self.Trkseg(
                        [self.Point(xy, {}, {}) for xy in track[0].vertices],
                        {"name":"segment0"}, {})], {}, {})
        self.assertEqual(g.tracks[0], expected)
        return

    def test_add_route(self):
        route = vector.guppy.Line([(np.random.random(), np.random.random())
                             for i in range(10)], properties={"name":"route0"})
        g = vector.gpx.GPX()
        g.add_route(route)
        expected = self.Route([self.Point(xy, {}, {}) for xy in route.vertices],
                              {"name":"route0"}, {})
        self.assertEqual(g.routes[0], expected)
        return




if __name__ == "__main__":
    unittest.main()

