import unittest
import numpy as np

import karta
import karta.vector.rtree
from karta.vector.rtree import RTree

class TestGeom(object):
    def __init__(self, bbox):
        self.bbox = bbox

class RTreeTests(unittest.TestCase):

    def test_build(self):
        np.random.seed(49)
        geoms = []
        for i in range(1000):
            offset = np.random.rand(1, 2)*5
            line = karta.Line(np.random.rand(5, 2)+offset)
            geoms.append(line)

        tree = RTree(geoms, maxchildren=25)
        self.assertEqual(type(tree), RTree)
        return

    def test_search_within(self):
        np.random.seed(49)
        geoms = []
        for i in range(1000):
            offset = np.random.rand(1, 2)*5
            line = karta.Line(np.random.rand(5, 2)+offset)
            geoms.append(line)

        rtree = RTree(geoms, maxchildren=25)
        corner_geoms = rtree.search_within((4, 4, 6, 6))
        self.assertEqual(len(corner_geoms), 50)
        return

    def test_search_overlapping(self):
        np.random.seed(49)
        geoms = []
        for i in range(1000):
            offset = np.random.rand(1, 2)*5
            line = karta.Line(np.random.rand(5, 2)+offset)
            geoms.append(line)

        rtree = RTree(geoms, maxchildren=25)
        corner_geoms = rtree.search_overlapping((4, 4, 6, 6))
        self.assertEqual(len(corner_geoms), 143)
        return

if __name__ == "__main__":
    unittest.main()
