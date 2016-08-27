import unittest
from karta import Point, Polygon
from karta.crs import LonLatWGS84, SphericalEarth
from karta.vector.dateline import crosses_dateline

class TestDateline(unittest.TestCase):

    def test_crosses_dateline(self):
        self.assertTrue(crosses_dateline(179, -179))
        self.assertTrue(crosses_dateline(-150, 140))
        self.assertFalse(crosses_dateline(160, 10))
        self.assertFalse(crosses_dateline(-20, 30))
        return

    def test_azimuth(self):
        for crs in (SphericalEarth, LonLatWGS84):
            pt0 = Point((0.0, 0.0), crs=crs)
            pt1 = Point((-1.0, 1.0), crs=crs)

            pt2 = Point((-179.5, 0.0), crs=crs)
            pt3 = Point((179.5, 1.0), crs=crs)
            self.assertAlmostEqual(pt0.azimuth(pt1), pt2.azimuth(pt3), places=8)

    def test_distance(self):
        for crs in (SphericalEarth, LonLatWGS84):
            pt0 = Point((0.0, 0.0), crs=crs)
            pt1 = Point((-1.0, 1.0), crs=crs)

            pt2 = Point((-179.5, 0.0), crs=crs)
            pt3 = Point((179.5, 1.0), crs=crs)
            self.assertAlmostEqual(pt0.distance(pt1), pt2.distance(pt3), places=8)

    def test_area(self):
        for crs in (SphericalEarth, LonLatWGS84):
            poly0 = Polygon([(-1, -1), (1, -1), (1, 1), (-1, 1)], crs=crs)
            poly1 = Polygon([(179, -1), (-179, -1), (-179, 1), (179, 1)], crs=crs)
            self.assertAlmostEqual(poly0.area, poly1.area)

    def test_bbox_geographical(self):
        for crs in (SphericalEarth, LonLatWGS84):
            poly = Polygon([(179, -1), (-179, -1), (-179, 1), (179, 1)], crs=crs)
            bb = poly.bbox
            self.assertEqual((bb[0], bb[2]), (179, -179))
            self.assertAlmostEqual(bb[1], -1.000152297, places=8)
            self.assertAlmostEqual(bb[3], 1.000152297, places=8)

if __name__ == "__main__":
    unittest.main()
