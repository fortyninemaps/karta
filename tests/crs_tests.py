
import unittest
import math
import karta.crs as crs

class TestCRS(unittest.TestCase):

    def test_CartesianProj(self):
        xp, yp = crs.Cartesian.proj(3, 4)
        self.assertEqual((xp, yp), (3, 4))
        return

    def test_CartesianGeodFwd(self):
        az = math.atan(0.75) * 180 / math.pi
        lons, lats, backaz = crs.Cartesian.geod.fwd(0.0, 0.0, az, 5)
        self.assertAlmostEqual(lons, 4.0, places=12)
        self.assertAlmostEqual(lats, 3.0, places=12)
        self.assertAlmostEqual(backaz, az+180.0, places=12)
        return

    def test_CartesianGeodInv(self):
        lon0 = 367
        lat0 = 78
        lon1 = 732
        lat1 = 23
        a1 = lon1-lon0
        a2 = lat1-lat0
        d_ = math.sqrt(a1**2 + a2**2)
        az_ = math.atan(a2/a1) + 2*math.pi
        baz_ = az_ - math.pi
        az, baz, d = crs.Cartesian.geod.inv(lon0, lat0, lon1, lat1, radians=True)
        self.assertAlmostEqual(d_, d, places=12)
        self.assertAlmostEqual(az_, az, places=12)
        self.assertAlmostEqual(baz_, baz, places=12)
        return

    def test_equal1(self):
        WGS84 = crs.CustomCRS(proj="+proj=longlat +datum=WGS84 +no_defs",
                              geod="+ellps=WGS84", name="WGS84 (Geographical)")
        WGS84_ = crs.CustomCRS(proj="+proj=longlat +datum=WGS84 +no_defs",
                               geod="+ellps=WGS84", name="WGS84 (Geographical)")
        self.assertTrue(WGS84 == WGS84_)
        return

    def test_equal2(self):
        WGS84 = crs.CustomCRS(proj="+proj=longlat +datum=WGS84 +no_defs",
                              geod="+ellps=WGS84", name="WGS84 (Geographical)")
        NAD83 = crs.CustomCRS(proj="+proj=longlat +datum=WGS84 +no_defs",
                              geod="+ellps=GRS80", name="NAD83 (Geographical)")
        self.assertTrue(not WGS84 == NAD83)
        return

    def test_not_equal1(self):
        WGS84 = crs.CustomCRS(proj="+proj=longlat +datum=WGS84 +no_defs",
                              geod="+ellps=WGS84", name="WGS84 (Geographical)")
        NAD83 = crs.CustomCRS(proj="+proj=longlat +datum=WGS84 +no_defs",
                              geod="+ellps=GRS80", name="NAD83 (Geographical)")
        self.assertTrue(WGS84 != NAD83)
        return

    def test_not_equal2(self):
        WGS84 = crs.CustomCRS(proj="+proj=longlat +datum=WGS84 +no_defs",
                              geod="+ellps=WGS84", name="WGS84 (Geographical)")
        WGS84_ = crs.CustomCRS(proj="+proj=longlat +datum=WGS84 +no_defs",
                               geod="+ellps=WGS84", name="WGS84 (Geographical)")
        self.assertTrue(not WGS84 != WGS84_)
        return

if __name__ == "__main__":
    unittest.main()

