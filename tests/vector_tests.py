""" Unit tests for vector functions """

import unittest
import os
import numpy as np
import karta.vector as vector
from karta.vector.geojson import GeoJSONReader
from test_helper import md5sum

class TestGuppy(unittest.TestCase):

    def setUp(self):
        self.point = vector.guppy.Point((1.0, 2.0, 3.0))

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
        return

    def test_point_bearing(self):
        other = vector.guppy.Point((7.0, 8.0))
        self.assertEqual(self.point.bearing(other), 225.0 / 180.0 * np.pi)
        return

    def test_point_azimuth(self):
        other = vector.guppy.Point((2.0, 1.0, 1.0))
        self.assertEqual(self.point.azimuth(other), -np.arctan(2./np.sqrt(2.0)))
        return

    def test_point_shift(self):
        point = vector.guppy.Point((-3.0, 5.0, 2.5))
        point.shift((4.0, -3.0, 0.5))
        self.assertEqual(self.point, point)
        return

    def test_multipoint_getset(self):
        self.assertEqual(self.mp[0], self.vertices[0])
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
        with open('data/testmp2vtp.vtp', 'w') as f:
            self.mp.to_vtk(f)
        self.assertEqual(md5sum('data/testmp2vtp.vtp'),
                         md5sum('reference_data/testmp2vtp.vtp'))
        return

    def test_geojson(self):
        # Test GeoJSON output for a Multipoint
        with open('data/testgeojson.json', 'w') as f:
            self.mp.to_geojson(f)
        self.assertEqual(md5sum('data/testgeojson.json'),
                         md5sum('reference_data/testgeojson.json'))
        return

class TestGeoJSONInput(unittest.TestCase):

    def test_point_read(self):
        with open('reference_data/point.geojson') as f:
            reader = GeoJSONReader(f)
        res = reader.pull_points()
        self.assertEqual(res[0].coordinates, [100.0, 0.0])
        return

    def test_linestring_read(self):
        with open('reference_data/linestring.geojson') as f:
            reader = GeoJSONReader(f)
        res = reader.pull_lines()
        self.assertEqual(res[0].coordinates, [[100.0, 0.0], [101.0, 1.0]])
        return

    def test_polygon_read(self):
        with open('reference_data/polygon.geojson') as f:
            reader = GeoJSONReader(f)
        res = reader.pull_polygons()
        self.assertEqual(res[0].coordinates, [[[100.0, 0.0], [101.0, 0.0],
                                               [101.0, 1.0], [100.0, 1.0],
                                               [100.0, 0.0]]])
        return

    def test_multipoint_read(self):
        with open('reference_data/multipoint.geojson') as f:
            reader = GeoJSONReader(f)
        res = reader.pull_multipoints()
        self.assertEqual(res[0].coordinates, [[100.0, 0.0], [101.0, 1.0]])
        return

    def test_multilinestring_read(self):
        with open('reference_data/multilinestring.geojson') as f:
            reader = GeoJSONReader(f)
        res = reader.pull_lines()
        self.assertEqual(res[0].coordinates, [[100.0, 0.0], [101.0, 1.0]])
        self.assertEqual(res[1].coordinates, [[102.0, 2.0], [103.0, 3.0]])
        return

    # This test will fail until holes are implemented
    def test_multipolygon_read(self):
        with open('reference_data/multipolygon.geojson') as f:
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
        self.points = [vector.gpx.Point((np.random.random(), np.random.random()), {}, {}) for i in range(20)]
        self.segments = [vector.gpx.Trkseg(self.points, {}, {})]
        self.tracks = [vector.gpx.Track(self.segments, {}, {})]
        self.routes = [vector.gpx.Route(self.points, {}, {})]
        return

    def test_track_init(self):
        """ Test initiation and writing of a GPX file containing a single track. """
        #g = vector.gpx.GPX(waypoints=self.points, tracks=self.tracks, routes=self.routes)
        g = vector.gpx.GPX()
        for i, pt in enumerate(self.points):
            g.waypts[i] = pt
        g.tracks = {0:self.tracks[0]}
        g.routes = {0:self.routes[0]}
        g.writefile("test.gpx")
        return

if __name__ == "__main__":
    unittest.main()

