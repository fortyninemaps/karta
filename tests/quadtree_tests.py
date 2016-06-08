""" Unit tests for vector functions """

import unittest
import numpy as np

import karta
from karta.vector.coordstring import CoordString
from karta.vector.quadtree import QuadTree
#from karta.vector.quadtree import QuadTree, Node
#from karta.vector.quadtree import addpt, split, hashpt
#from karta.vector.quadtree import querypt_recursion, querypt_hash
#from karta.vector.quadtree import iswithin, overlaps

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
