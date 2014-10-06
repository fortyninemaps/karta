""" Test constructing guppy geometry instances """

import unittest
import numpy as np
import karta
from karta import Point, Multipoint, Line, Polygon

LONLAT = karta.crs.LONLAT

class TestGuppy(unittest.TestCase):

    def setUp(self):
        x = range(-5, 5)
        y = [x_**2 for x_ in x]
        self.vertices = list(zip(x, y))
        self.data = [x_*y_ for (x_, y_) in self.vertices]
        return

    def test_multipoint_from_points(self):
        pts = [Point((x, y), data={"d": d}, properties={"p":i}, crs=LONLAT)
                for i,((x,y),d) in enumerate(zip(self.vertices, self.data))]

        mp = Multipoint(pts)
        ans = Multipoint(self.vertices, data={"d":self.data}, crs=LONLAT)
        self.assertEqual(mp, ans)
        return

    def test_line_from_points(self):
        pts = [Point((x, y), data={"d": d}, properties={"p":i}, crs=LONLAT)
                for i,((x,y),d) in enumerate(zip(self.vertices, self.data))]

        mp = Line(pts)
        ans = Line(self.vertices, data={"d":self.data}, crs=LONLAT)
        self.assertEqual(mp, ans)
        return

    def test_polygon_from_points(self):
        pts = [Point((x, y), data={"d": d}, properties={"p":i}, crs=LONLAT)
                for i,((x,y),d) in enumerate(zip(self.vertices, self.data))]

        mp = Polygon(pts)
        ans = Polygon(self.vertices, data={"d":self.data}, crs=LONLAT)
        self.assertEqual(mp, ans)
        return

if __name__ == "__main__":
    unittest.main()

