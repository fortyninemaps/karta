""" Test constructing geometry instances """
import unittest
import numpy as np
from karta import Point, Line, Polygon, Multipoint, Multiline, Multipolygon
from karta.crs import LonLatWGS84, WebMercator

class MergeTests(unittest.TestCase):

    def test_merge_points(self):
        pt1 = Point((1, 2), properties={"number": 16})
        mp1 = Multipoint([(5, 6), (5, 4), (4, 3)],
                         data={"number": [2, 76, 4]})
        pt2 = Point((3, 4, 5), properties={"number": 22, "temperature": "hot"})
        mp2 = Multipoint([(5, 6), (4, 7), (3, 9), (4, 9)],
                         data={"number": [4, 674, 32, 56],
                               "temperature": ["cool", "hot", "scorching", "mild"]})

        merged = Multipoint.merge(pt1, mp1, pt2, mp2)
        self.assertEqual(len(merged), 9)
        self.assertEqual(merged.vertices.asarray().shape, (9, 2))
        self.assertEqual(merged.d["number"], [16, 2, 76, 4, 22, 4, 674, 32, 56])

    def test_merge_lines(self):
        ml1 = Multiline([[(1, 2), (3, 4), (5, 6)], [(3, 4), (2, 3)]])
        ln1 = Line([(5, 6, 1), (5, 4, 2), (4, 3, 0)])
        ml2 = Multiline([[(4, 2), (8, 4), (1, 5)], [(2, 4), (8, 6)]])
        ln2 = Line([(5, 6), (4, 7), (3, 9), (4, 9)])

        merged = Multiline.merge(ml1, ln1, ml2, ln2)
        self.assertEqual(len(merged), 6)
        self.assertEqual(merged.vertices[0].asarray().shape, (3, 2))
        self.assertEqual(merged.vertices[1].asarray().shape, (2, 2))

    def test_merge_polygons(self):
        pg1 = Polygon([[(3, 4), (6, 4), (7, 6), (2, 5)]])
        pg2 = Polygon([[(2, 4), (5, 4), (6, 6), (0, 5)]])
        mp1 = Multipolygon([[[(0, 0), (9, 1), (5, 8)], [(1, 1), (6, 2), (3, 3)]], [[(3, 3), (4, 3), (4, 5), (4, 3)]]])
        mp2 = Multipolygon([[[(1, 9), (6, 0), (9, 8)], [(4, 3), (4, 4), (5, 4)]], [[(6, 2), (7, 3), (8, 4), (7, 5)]]])

        merged = Multipolygon.merge(pg1, mp1, pg2, mp2)
        self.assertEqual(len(merged), 6)

if __name__ == "__main__":
    unittest.main()

