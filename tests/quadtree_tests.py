""" Unit tests for vector functions """

import unittest
import ctypes
from ctypes import c_int, c_double, c_char_p, c_bool
import numpy as np

import karta
from karta.vector.coordstring import CoordString
from karta.vector.quadtree import QuadTree

class BBOX(ctypes.Structure):
    _fields_ = [("xmin", c_double), ("ymin", c_double),
                ("xmax", c_double), ("ymax", c_double)]

class POSITION(ctypes.Structure):
    _fields_ = [("id", c_int), ("x", c_double), ("y", c_double)]

class TestCQuadTree(unittest.TestCase):

    def setUp(self):
        self.qt = ctypes.CDLL(karta.vector.quadtree.__file__)

    def test_mind(self):
        self.qt.mind.restype = c_double
        self.assertEqual(self.qt.mind(c_double(3.0), c_double(4.0)), 3.0)
        self.assertEqual(self.qt.mind(c_double(3.0), c_double(-2.0)), -2.0)
        return

    def test_maxd(self):
        self.qt.maxd.restype = c_double
        self.assertEqual(self.qt.maxd(c_double(3.0), c_double(4.0)), 4.0)
        self.assertEqual(self.qt.maxd(c_double(-3.0), c_double(-4.0)), -3.0)
        return

    def test_iswithin(self):
        self.qt.iswithin.restype = c_bool
        self.qt.qt_new_bbox.restype = ctypes.POINTER(BBOX)
        self.qt.qt_new_position.restype = ctypes.POINTER(POSITION)
        bboxPtr = self.qt.qt_new_bbox(c_double(-1.0), c_double(-1.0),
                                      c_double(1.0), c_double(1.0))
        ptPtr = self.qt.qt_new_position(c_int(0), c_double(0.5), c_double(0.5))
        self.assertTrue(self.qt.iswithin(bboxPtr, ptPtr))

        ptPtr = self.qt.qt_new_position(c_int(1), c_double(-1.5), c_double(0.5))
        self.assertFalse(self.qt.iswithin(bboxPtr, ptPtr))

        ptPtr = self.qt.qt_new_position(c_int(1), c_double(0.5), c_double(1.5))
        self.assertFalse(self.qt.iswithin(bboxPtr, ptPtr))
        return

class TestQuadTree(unittest.TestCase):

    def test_quadtree_construction(self):
        vertices = [(x, x**2) for x in np.linspace(-10, 10, 1000)]
        cs = CoordString(vertices)
        quadtree = QuadTree(cs)
        self.assertEqual(len(quadtree), 1000)
        return

    def test_quadtree_search(self):
        vertices = [(x, x**2) for x in np.linspace(-10, 10, 1000)]
        cs = CoordString(vertices)
        quadtree = QuadTree(cs)

        indices_within = quadtree.search_within(-5, 0, 5, 25)
        self.assertEqual(len(indices_within), 500)
        self.assertEqual(min(indices_within), 250)
        self.assertEqual(max(indices_within), 749)
        return

    def test_quadtree_duplicates(self):
        # a naive quadtree with enter an infinite loop when there are
        # duplicates exceeding leaf node capacity. this test ensures that this
        # doesn't happen, and that the duplicate points can be retrieved
        vertices = [(3.0, 4.0) for _ in range(11)]
        cs = CoordString(vertices)
        quadtree = QuadTree(cs, leaf_capacity=10)
        indices_within = quadtree.search_within(2, 3, 4, 5)
        self.assertEqual(len(indices_within), 11)
        return

class TestGeometryWithQuadTree(unittest.TestCase):

    def test_within_radius(self):
        """ Get the points near to another point using a quadtree approach. """
        crs = karta.crs.LonLatWGS84
        pt = karta.Point((0.0, 30.0), crs=crs)
        np.random.seed(42)
        x = (np.random.random(1000) - 0.5) * 180.0
        y = (np.random.random(1000) - 0.5) * 30.0
        mp_noindex = karta.Multipoint(zip(x, y), crs=crs, build_index=False)
        mp = karta.Multipoint(zip(x, y), crs=crs, build_index=True)

        subset_noindex = mp_noindex.within_radius(pt, 5e6)
        subset = mp.within_radius(pt, 5e6)
        self.assertEqual(len(subset), len(subset_noindex))
        for i in range(0, len(subset), 50):
            self.assertTrue(subset[i] in subset_noindex)
        return

    def test_polygon_contains(self):
        """ Search for points contained within a polygon, using a quadtree. """
        crs = karta.crs.SphericalEarth
        th = np.linspace(-np.pi, np.pi, 18)
        xp = 5*np.cos(th)
        yp = 5*np.sin(th)
        poly = karta.Polygon(zip(xp, yp), crs=crs)

        np.random.seed(42)
        x = (np.random.random(1000) - 0.5) * 180.0
        y = (np.random.random(1000) - 0.5) * 30.0
        mp_noindex = karta.Multipoint(zip(x, y), crs=crs, build_index=False)
        mp = karta.Multipoint(zip(x, y), crs=crs, build_index=True)

        contained_noindex = mp_noindex.within_polygon(poly)
        contained = mp.within_polygon(poly)
        self.assertEqual(len(contained), len(contained_noindex))
        for i in range(len(contained)):
            self.assertTrue(contained[i] in contained_noindex)
        return

if __name__ == "__main__":
    unittest.main()
