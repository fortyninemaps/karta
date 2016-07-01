import unittest
import ctypes
from ctypes import c_int, c_float, pointer
import numpy as np

import karta
import karta.vector.rtree
from karta.vector.rtree import RTree

class TestGeom(object):
    def __init__(self, bbox):
        self.bbox = bbox

class BBOX(ctypes.Structure):
    _fields_ = [("xmin", c_float), ("ymin", c_float),
                ("xmax", c_float), ("ymax", c_float)]

class CRTreeTests(unittest.TestCase):

    def setUp(self):
        self.rt = ctypes.CDLL(karta.vector.rtree.__file__)

    def test_is_within(self):
        self.rt.is_within.restype = c_int
        bbox1 = pointer(BBOX(0.0, 0.0, 2.0, 2.0))
        bbox2 = pointer(BBOX(0.5, 0.5, 1.5, 1.5))
        self.assertEqual(self.rt.is_within(bbox1, bbox2), 1)

        bbox1 = pointer(BBOX(0.6, 0.0, 2.0, 2.0))
        bbox2 = pointer(BBOX(0.5, 0.5, 1.5, 1.5))
        self.assertEqual(self.rt.is_within(bbox1, bbox2), 0)

        bbox1 = pointer(BBOX(0.0, 0.0, 2.0, 2.0))
        bbox2 = pointer(BBOX(0.1, 0.5, 1.1, 1.5))
        self.assertEqual(self.rt.is_within(bbox1, bbox2), 1)
        return

    def test_is_overlapping(self):
        self.rt.is_overlapping.restype = c_int
        bbox1 = pointer(BBOX(0.0, 0.0, 1.0, 1.0))
        bbox2 = pointer(BBOX(0.5, 0.5, 1.5, 1.5))
        self.assertEqual(self.rt.is_overlapping(bbox1, bbox2), 1)

        bbox1 = pointer(BBOX(0.0, 0.0, 1.0, 1.0))
        bbox2 = pointer(BBOX(1.5, 1.5, 2.5, 2.5))
        self.assertEqual(self.rt.is_overlapping(bbox1, bbox2), 0)

        bbox1 = pointer(BBOX(0.0, 0.0, 2.0, 2.0))
        bbox2 = pointer(BBOX(-0.5, 0.0, 0.5, 0.8))
        self.assertEqual(self.rt.is_overlapping(bbox1, bbox2), 1)
        return

    def test_bbox_intersection_area(self):
        self.rt.bbox_intersection_area.restype = c_float
        bbox1 = pointer(BBOX(0.0, 0.0, 1.0, 1.0))
        bbox2 = pointer(BBOX(0.5, 0.5, 1.5, 1.5))
        self.assertAlmostEqual(self.rt.bbox_intersection_area(bbox1, bbox2), 0.25, places=6)

        bbox1 = pointer(BBOX(0.0, 0.0, 1.0, 1.0))
        bbox2 = pointer(BBOX(1.5, 1.5, 2.5, 2.5))
        self.assertAlmostEqual(self.rt.bbox_intersection_area(bbox1, bbox2), 0.0, places=6)

        bbox1 = pointer(BBOX(0.0, 0.0, 2.0, 2.0))
        bbox2 = pointer(BBOX(-0.5, 0.0, 0.5, 0.8))
        self.assertAlmostEqual(self.rt.bbox_intersection_area(bbox1, bbox2), 0.4, places=6)
        return

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
