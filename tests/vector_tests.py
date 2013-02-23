""" Unit tests for vector functions """

import unittest
import os
from test_helper import md5sum
import geo_tools.vector as vector

class TestGuppy(unittest.TestCase):

    def setUp(self):
        self.poly = vector.guppy.Polygon([(0.0, 8.0), (0.0, 5.0), (6.0, 1.0)])


    def test_point_creation(self):
        # Create a point
        P = vector.guppy.Point((1.0, 2.0, 3.0))
        return

    def test_multipoint_create(self):
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

        L0 = vector.guppy.Multipoint(vertices)
        L1 = vector.guppy.Multipoint(vertices, data=data0)
        L2 = vector.guppy.Multipoint(vertices, data={'d0':data0, 'd1':data1})
        return

    def test_poly_vertices(self):
        self.assertEqual(self.poly.get_vertices(),
                         [(0.0, 8.0), (0.0, 5.0), (6.0, 1.0)])
        return

    def test_poly_coordinates(self):
        self.assertEqual(self.poly.get_coordinate_lists(),
                         ([0.0, 0.0, 6.0], [8.0, 5.0, 1.0], [0.0, 0.0, 0.0]))
        return

    def test_poly_extents(self):
        self.assertEqual(self.poly.get_extents(), [(0.0, 6.0), (1.0, 8.0)])
        return

    def test_poly_length(self):
        self.assertEqual(self.poly.length(), 10.21110255092798)
        return


class TestGuppyIO(unittest.TestCase):

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


if __name__ == "__main__":
    unittest.main()

