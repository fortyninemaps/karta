
import unittest
import math
from karta.crs import crsreg
import karta.crs2 as crs2


class TestCRS(unittest.TestCase):

    def test_CartesianProj(self):
        xp, yp = crs2.Cartesian.proj(3, 4)
        self.assertEqual((xp, yp), (3, 4))
        return

    def test_CartesianGeodFwd(self):
        az = math.atan(0.75) * 180 / math.pi
        lons, lats, backaz = crs2.Cartesian.geod.fwd(0.0, 0.0, az, 5)
        self.assertEqual(lons, 4.0)
        self.assertEqual(lats, 3.0)
        self.assertEqual(backaz, az+180.0)
        return

    def test_CartesianGeodInv(self):
        raise NotImplementedError()
        return

    def test_equal1(self):
        WGS84 = crsreg.get_CRS("LONLAT_WGS84")
        WGS84_ = crsreg.get_CRS("LONLAT_WGS84")
        self.assertTrue(WGS84 == WGS84_)
        return

    def test_equal2(self):
        WGS84 = crsreg.get_CRS("LONLAT_WGS84")
        NAD83 = crsreg.get_CRS("LONLAT_NAD83")
        self.assertTrue(not WGS84 == NAD83)
        return

    def test_not_equal1(self):
        WGS84 = crsreg.get_CRS("LONLAT_WGS84")
        NAD83 = crsreg.get_CRS("LONLAT_NAD83")
        self.assertTrue(WGS84 != NAD83)
        return

    def test_not_equal2(self):
        WGS84 = crsreg.get_CRS("LONLAT_WGS84")
        WGS84_ = crsreg.get_CRS("LONLAT_WGS84")
        self.assertTrue(not WGS84 != WGS84_)
        return

if __name__ == "__main__":
    unittest.main()

