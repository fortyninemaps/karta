""" Test constructing geometry instances """

import unittest
from karta import Point, Multipoint, Line, Polygon
from karta.vector.geometry import multipart_from_singleparts
from karta.crs import LonLatWGS84

class TestVectorGeometry(unittest.TestCase):

    def setUp(self):
        x = range(-5, 5)
        y = [x_**2 for x_ in x]
        self.vertices = list(zip(x, y))
        self.data = [x_*y_ for (x_, y_) in self.vertices]
        return

    def test_multipart_from_points(self):
        pts = [Point((x, y), properties={"p":i}, crs=LonLatWGS84)
                for i,((x,y),d) in enumerate(zip(self.vertices, self.data))]

        mp = multipart_from_singleparts(pts)
        ans = Multipoint(self.vertices, data={"p":range(len(pts))}, crs=LonLatWGS84)
        self.assertEqual(mp, ans)
        return

    def test_line_from_points(self):
        pts = [Point((x, y), properties={"p":i}, crs=LonLatWGS84)
                for i,((x,y),d) in enumerate(zip(self.vertices, self.data))]

        line = Line([pt.vertex for pt in pts], crs=LonLatWGS84)
        ans = Line(self.vertices, crs=LonLatWGS84)
        self.assertEqual(line, ans)
        return

    def test_polygon_from_points(self):
        pts = [Point((x, y), properties={"p":i}, crs=LonLatWGS84)
                for i,((x,y),d) in enumerate(zip(self.vertices, self.data))]
        poly = Polygon([pt.vertex for pt in pts], crs=LonLatWGS84)
        ans = Polygon(self.vertices, crs=LonLatWGS84)
        self.assertEqual(poly, ans)
        return

    def test_multipoint_datadict(self):
        # create a line
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

        Multipoint(vertices, data={'d0':data0, 'd1':data1})
        return

if __name__ == "__main__":
    unittest.main()
