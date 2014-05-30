""" Unit tests for vector functions """

import unittest
import numpy as np

import karta
from karta.vector.quadtree import QuadTree, Node
from karta.vector.quadtree import addpt, split, hashpt
from karta.vector.quadtree import querypt_recursion, querypt_hash
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
        (node, _) = addpt(node, (2.5, 2.5), 1, 25, 100)
        self.assertEqual(node.children[0].children, 
                         [(x,y) for x in range(3) for y in range(3)])
        self.assertEqual(node.children[1].children,
                         [(x,y) for x in range(3,5) for y in range(3)])
        self.assertEqual(node.children[2].children,
                         [(x,y) for x in range(3) for y in range(3,5)])
        self.assertEqual(node.children[3].children,
                         [(x,y) for x in range(3,5) for y in range(3,5)] + [(2.5,2.5)])

    def test_hashpt1(self):
        # Find the depth at which (0.99, 0.99) is not in the last quadrant (should be 7)
        hashgen = hashpt(0, 1, 0, 1, 0.99, 0.99)
        for (i,h) in enumerate(hashgen):
            if h == 0:
                break
            if i == 10:
                break
        self.assertEqual(i, 6)
        return

    def test_hashpt2(self):
        # Find the hash sequence for (0.26, 0.84) out to length 10
        hashgen = hashpt(0, 1, 0, 1, 0.26, 0.84)
        hsh = [next(hashgen) for i in range(10)]
        self.assertEqual(hsh, [2, 3, 0, 2, 0, 2, 3, 2, 1, 0])
        return

    def test_querypt_recursion(self):
        pts = [(x**0.5,0.5*y**0.875) for x in range(50) for y in range(50)]
        node = Node([], (0, 7.5, 0, 16), True)
        for pt in pts:
            node,_ = addpt(node, pt, 1, 20, 200)
        testpts = [(x**0.5,0.5*y**0.875) for (x,y) in zip((3,12,44,23,36),
                                                          (46,42,28,2,13))]
        shouldbetrue = [querypt_recursion(node, testpt) for testpt in testpts]
        testpts = [(x**0.5,0.5*y**0.875) for (x,y) in zip((73,12,54,23,63),
                                                          (46,72,28,82,13))]
        shouldbefalse = [querypt_recursion(node, testpt) for testpt in testpts]
        self.assertTrue(False not in shouldbetrue)
        self.assertTrue(True not in shouldbefalse)
        return

    def test_querypt_hash(self):
        pts = [(x**0.5,0.5*y**0.875) for x in range(50) for y in range(50)]
        node = Node([], (0, 7.5, 0, 16), True)
        for pt in pts:
            node,_ = addpt(node, pt, 1, 20, 200)
        testpts = [(x**0.5,0.5*y**0.875) for (x,y) in zip((3,12,44,23,36),
                                                          (46,42,28,2,13))]
        shouldbetrue = [querypt_hash(node, testpt) for testpt in testpts]
        testpts = [(x**0.5,0.5*y**0.875) for (x,y) in zip((73,12,54,23,63),
                                                          (46,72,28,82,13))]
        shouldbefalse = [querypt_hash(node, testpt) for testpt in testpts]
        self.assertTrue(False not in shouldbetrue)
        self.assertTrue(True not in shouldbefalse)
        return

    def test_querypt_big(self):
        pts = [(x**0.5,0.05*y**0.875) for x in range(500) for y in range(500)]
        node = Node([], (0, 32, 0, 22), True)
        for pt in pts:
            node,_ = addpt(node, pt, 1, 20, 200)
        testpts = [(x**0.5,0.05*y**0.875) for (x,y) in zip((3,12,44,23,36),
                                                          (46,42,28,2,13))]
        shouldbetrue = [querypt_recursion(node, testpt) for testpt in testpts]
        testpts = [(x**0.5,0.05*y**0.875) for (x,y) in zip((730,128,540,231,630),
                                                          (46,720,28,820,13))]
        shouldbefalse = [querypt_recursion(node, testpt) for testpt in testpts]
        self.assertTrue(False not in shouldbetrue)
        self.assertTrue(True not in shouldbefalse)
        return

    def test_quadtree_add_query(self):
        pts = [(x**0.5,0.5*y**0.875) for x in range(50) for y in range(50)]
        QTree = QuadTree((0, 32, 0, 22))
        for pt in pts:
            QTree.addpt(pt)
        testpts = [(x**0.5,0.5*y**0.875) for (x,y) in zip((3,12,44,23,36),
                                                          (46,42,28,2,13))]
        shouldbetrue = [QTree.querypt(pt) for pt in testpts]
        testpts = [(x**0.5,0.5*y**0.875) for (x,y) in zip((73,12,54,23,63),
                                                          (46,72,28,82,13))]
        shouldbefalse = [QTree.querypt(pt) for pt in testpts]
        self.assertTrue(False not in shouldbetrue)
        self.assertTrue(True not in shouldbefalse)

    def test_quadtree_add_bbox(self):
        pts = [(x**0.5,0.5*y**0.875) for x in range(50) for y in range(50)]
        QTree = QuadTree((0, 32, 0, 22))
        for pt in pts:
            QTree.addpt(pt)
        pts = QTree.getfrombbox((16, 24, 8, 12))
        ans = [pt for pt in pts if (16 <= pt[0] < 24) and (8 <= pt[1] < 12)]
        pts.sort()
        ans.sort()
        self.assertEqual(pts, ans)
        return


if __name__ == "__main__":
    unittest.main()

