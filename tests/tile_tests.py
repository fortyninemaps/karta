import unittest
from karta.vector.geometry import get_tile_tuple, tile_nw_corner, tile_bbox
from karta.vector.geometry import Point
from karta.crs import LonLatWGS84

class TileTests(unittest.TestCase):

    def test_get_tile_tuple(self):
        t = get_tile_tuple(Point((0, 0), crs=LonLatWGS84), 0)
        self.assertEqual(t, (0, 0, 0))
        
        t = get_tile_tuple(Point((0, 0), crs=LonLatWGS84), 8)
        self.assertEqual(t, (8, 128, 128))

        t = get_tile_tuple(Point((60, -30), crs=LonLatWGS84), 12)
        self.assertEqual(t, (12, 2730, 2406))
        return

    def test_tile_nw_corner(self):
        pt = tile_nw_corner(0, 0, 0)
        self.assertEqual(pt.x, -180.0)
        self.assertAlmostEqual(pt.y, 85.05112877, places=7)

        pt = tile_nw_corner(1, 1, 1)
        self.assertEqual(pt.vertex, (0, 0))
        return

    def test_tile_bbox(self):
        bbox = tile_bbox(1, 0, 0)
        self.assertEqual(bbox[0], -180.0)
        self.assertEqual(bbox[2], 0.0)
        self.assertEqual(bbox[1], 0.0)
        self.assertAlmostEqual(bbox[3], 85.05112877, places=7)
        return

if __name__ == "__main__":
    unittest.main()