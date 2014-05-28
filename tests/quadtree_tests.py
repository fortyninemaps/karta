""" Unit tests for vector functions """

import unittest
import numpy as np

import karta
from karta.vector.quadtree import QuadTree, Node, addpt, split
from karta.vector.quadtree import iswithin, overlaps

class TestQuadTree(unittest.TestCase):

    def test_iswithin(self):
        bbox = (0, 4, 0, 4)
        self.assertTrue(iswithin(bbox, (2,2)))
        self.assertTrue(iswithin(bbox, (0,0)))
        self.assertFalse(iswithin(bbox, (5,6)))
        self.assertFalse(iswithin(bbox, (4,4)))
        return

    def test_addpt_leaf(self):
        node = Node([], (0, 1, 0, 1), True)
        node,_ = addpt(node, (0.5, 0.5), 1, 20, 100)
        self.assertEqual(node.children, [(0.5, 0.5)])
        return

    def test_split(self):
        pts = [(x,y) for x in range(5) for y in range(5)]
        node = Node(pts, (0, 4.1, 0, 4.1), True)
        parent = split(node)
        self.assertEqual(parent.children[0].children, 
                         [(x,y) for x in range(3) for y in range(3)])
        self.assertEqual(parent.children[1].children,
                         [(x,y) for x in range(3,5) for y in range(3)])
        self.assertEqual(parent.children[2].children,
                         [(x,y) for x in range(3) for y in range(3,5)])
        self.assertEqual(parent.children[3].children,
                         [(x,y) for x in range(3,5) for y in range(3,5)])
        return

    def test_split_bbox(self):
        pts = [(x,y) for x in range(5) for y in range(5)]
        node = Node(pts, (0, 4, 0, 4), True)
        parent = split(node)
        self.assertEqual(parent.children[0].bbox, (0, 2, 0, 2))
        self.assertEqual(parent.children[1].bbox, (2, 4, 0, 2))
        self.assertEqual(parent.children[2].bbox, (0, 2, 2, 4))
        self.assertEqual(parent.children[3].bbox, (2, 4, 2, 4))
        return

    def test_addpt_full(self):
        pts = [(x,y) for x in range(5) for y in range(5)]
        node = Node(pts, (0, 4.1, 0, 4.1), True)
        node,_ = addpt(node, (2.5, 2.5), 1, 25, 100)
        self.assertEqual(node.children[0].children, 
                         [(x,y) for x in range(3) for y in range(3)])
        self.assertEqual(node.children[1].children,
                         [(x,y) for x in range(3,5) for y in range(3)])
        self.assertEqual(node.children[2].children,
                         [(x,y) for x in range(3) for y in range(3,5)])
        self.assertEqual(node.children[3].children,
                         [(x,y) for x in range(3,5) for y in range(3,5)] + [(2.5,2.5)])



if __name__ == "__main__":
    unittest.main()

