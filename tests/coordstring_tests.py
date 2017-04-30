import unittest
import numpy as np
from karta.vector.coordstring import CoordString

class CoordStringTests(unittest.TestCase):

    def test_creation2(self):
        x = np.linspace(0, 5)
        y = x**2
        cs = CoordString(np.c_[x, y])
        self.assertEqual(cs.rank, 2)
        self.assertEqual(len(cs), len(x))
        return

    def test_creation3(self):
        x = np.linspace(0, 5)
        y = x**2
        z = x**3
        cs = CoordString(np.c_[x, y, z])
        self.assertEqual(cs.rank, 3)
        self.assertEqual(len(cs), len(x))
        return

    def test_creation_invalid(self):
        x = np.linspace(0, 5)
        y = x**2
        z = x**3
        with self.assertRaises(ValueError):
            # this one has too few dimensions
            cs = CoordString(x)
        with self.assertRaises(ValueError):
            # this one has too many dimensions
            cs = CoordString(np.c_[x, y, z, z])
        with self.assertRaises(ValueError):
            # this one has the wrong shape
            cs = CoordString(np.r_[x, y, z])
        return

    def test_bbox2(self):
        x = np.linspace(0, 5)
        y = x**2
        cs = CoordString(np.c_[x, y])
        self.assertEqual(cs.bbox, (0, 0, 5, 25))
        return

    def test_bbox3(self):
        x = np.linspace(0, 5)
        y = x**2
        z = x**3
        cs = CoordString(np.c_[x, y, z])
        self.assertEqual(cs.bbox, (0, 0, 5, 25))
        return

    def test_setitem(self):
        x = np.linspace(0, 5)
        y = x**2
        cs = CoordString(np.c_[x, y])
        cs[22] = np.array((-1, -3), dtype=np.float64)
        self.assertEqual(cs[22], (-1, -3))
        return

    def test_slicing(self):
        x = np.linspace(0, 5)
        y = x**2
        cs = CoordString(np.c_[x, y])
        self.assertTrue(np.all(cs.slice(10, 25) == np.c_[x, y][10:25]))
        return

    def test_slicing_stepped(self):
        x = np.linspace(0, 5)
        y = x**2
        cs = CoordString(np.c_[x, y])
        self.assertTrue(np.all(cs.slice(10, 35, 3) == np.c_[x, y][10:35:3]))
        return

    def test_asarray(self):
        x = np.linspace(0, 5)
        y = x**2
        cs = CoordString(np.c_[x, y])
        arr = cs.asarray()
        self.assertTrue(np.all(arr == np.c_[x, y]))
        return

    def test_asarray_empty(self):
        cs = CoordString([])
        arr = cs.asarray()
        self.assertTrue(np.all(arr == np.array([[]], dtype=np.float64)))

    def test_hash(self):
        A = CoordString([(i, i+1) for i in range(0, 100, 3)])
        B = CoordString([(i+1, i) for i in range(100, 0, -3)])
        C = CoordString([(i, i+1) for i in range(0, 100, 3)])
        self.assertEqual(hash(A), hash(C))
        self.assertNotEqual(hash(A), hash(B))

    def test_nan_raises(self):
        coords = [(1,2), (3,4), (5,np.nan), (7,8), (9,10)]
        with self.assertRaises(ValueError):
            CoordString(coords)

if __name__ == "__main__":
    unittest.main()
