import unittest
import numpy as np

import karta
import karta.vector.rtree as rtree

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

        tree = rtree.build_tree(geoms, max_children=25)
        self.assertEqual(type(tree), rtree.NonLeafNode)
        self.assertEqual(rtree.depth(tree), 3)
        self.assertEqual(len(tree.children), 4)
        return

    def test_leaf_split(self):
        node = rtree.LeafNode((0, 0, 1, 1), None, max_children=10)
        node.children = [TestGeom((0.1*i, 0.1*i, 0.1*i+0.05, 0.1*i+0.05))
                         for i in range(10)]
        newnodes = node.split()
        self.assertEqual(len(newnodes[0].children), 5)
        self.assertEqual(len(newnodes[1].children), 5)
        for child in node.children:
            self.assertTrue((child in newnodes[0].children) or
                            (child in newnodes[1].children))
        self.assertTrue(newnodes[0].parent is None)
        self.assertTrue(newnodes[1].parent is None)

    def test_nonleaf_split(self):
        node = rtree.NonLeafNode((0, 0, 1, 1), None, max_children=10)
        node.children = [rtree.LeafNode((0.1*i, 0.1*i, 0.1*i+0.05, 0.1*i+0.05), None)
                         for i in range(10)]
        newnodes = node.split()
        self.assertEqual(len(newnodes[0].children), 5)
        self.assertEqual(len(newnodes[1].children), 5)
        for child in node.children:
            self.assertTrue((child in newnodes[0].children) or
                            (child in newnodes[1].children))
        self.assertTrue(newnodes[0].parent is None)
        self.assertTrue(newnodes[1].parent is None)

if __name__ == "__main__":
    unittest.main()
