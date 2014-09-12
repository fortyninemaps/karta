
import unittest
from karta.crs import crsreg


class TestCRS(unittest.TestCase):

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

