
import unittest
import math
import numpy as np
import karta.crs as crs
import karta.geodesy as geodesy
import karta.errors

class TestCRS(unittest.TestCase):

    def test_get_proj4(self):
        c0 = crs.SphericalEarth
        self.assertEqual(c0.get_proj4(), "+proj=lonlat +a=6371009.000000 +b=6371009.000000 +no_defs")

        c1 = crs.LonLatWGS84
        proj4 = c1.get_proj4()
        self.assertTrue("+proj=lonlat" in proj4)
        self.assertTrue("+a=6378137.0" in proj4)
        self.assertTrue("+b=6356752.314245" in proj4)

        c2 = crs.NSIDCNorth
        proj4 = c2.get_proj4()
        self.assertTrue("+proj=stere" in proj4)
        self.assertTrue("+lat_0=90" in proj4)
        self.assertTrue("+lat_ts=70" in proj4)
        self.assertTrue("+lon_0=-45" in proj4)
        self.assertTrue("+k=1" in proj4)
        self.assertTrue("+x_0=0" in proj4)
        self.assertTrue("+y_0=0" in proj4)
        self.assertTrue("+units=m" in proj4)
        self.assertTrue("+datum=WGS84" in proj4)
        self.assertTrue("+a=6378137.0" in proj4)
        self.assertTrue("+f=0.0033528" in proj4)
        return

    def test_get_wkt(self):
        c0 = crs.SphericalEarth
        self.assertEqual(c0.get_wkt(), '+proj=lonlat +a=6371009.000000 +b=6371009.000000 +no_defs')

        c1 = crs.LonLatWGS84
        wkt = c1.get_wkt()
        self.assertTrue('GEOGCS["unnamed ellipse",DATUM["unknown",SPHEROID["unnamed",6378137,298.25722' in wkt)
        self.assertTrue('PRIMEM["Greenwich",0],UNIT["degree",0.01745329' in wkt)

        c2 = crs.NSIDCNorth
        wkt = c2.get_wkt()
        self.assertTrue('PROJCS["unnamed",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.25722' in wkt)
        self.assertTrue('AUTHORITY["EPSG","7030"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329' in wkt)
        self.assertTrue('AUTHORITY["EPSG","9108"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",70],PARAMETER["central_meridian",0],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]' in wkt)
        return

    def test_CartesianProj(self):
        xp, yp = crs.Cartesian.project(3.0, 4.0)
        self.assertEqual((xp, yp), (3.0, 4.0))
        return

    def test_CartesianForward(self):
        az = math.atan(0.75) * 180.0 / math.pi
        lons, lats, backaz = crs.Cartesian.forward(0.0, 0.0, az, 5)
        self.assertAlmostEqual(lons, 3.0, places=12)
        self.assertAlmostEqual(lats, 4.0, places=12)
        self.assertAlmostEqual(backaz, az+180.0, places=12)
        return

    def test_CartesianInverse(self):
        lon0 = 367.0
        lat0 = 78.0
        lon1 = 93.0
        lat1 = 23.0
        a1 = lon1-lon0
        a2 = lat1-lat0
        d_ = math.sqrt(a1**2 + a2**2)
        az_ = 1.5*math.pi - math.atan(a2/a1)
        baz_ = az_ - math.pi
        az, baz, d = crs.Cartesian.inverse(lon0, lat0, lon1, lat1, radians=True)
        self.assertAlmostEqual(d_, d, places=12)
        self.assertAlmostEqual(az_, az, places=12)
        self.assertAlmostEqual(baz_, baz, places=12)
        return

    def test_CartesianInverse2(self):
        lon0 = 367.0
        lat0 = 78.0
        lon1 = 732.0
        lat1 = 23.0
        a1 = lon1-lon0
        a2 = lat1-lat0
        d_ = math.sqrt(a1**2 + a2**2)
        az_ = 0.5*math.pi - math.atan(a2/a1)
        baz_ = az_ + math.pi
        az, baz, d = crs.Cartesian.inverse(lon0, lat0, lon1, lat1, radians=True)
        self.assertAlmostEqual(d_, d, places=12)
        self.assertAlmostEqual(az_, az, places=12)
        self.assertAlmostEqual(baz_, baz, places=12)
        return

    def test_SphericalForward1(self):
        lon1 = 0.0
        lat1 = 0.0
        lon2, lat2, baz = crs.SphericalEarth.forward(lon1, lat1, 90.0, 5003778.767588614,
                                                     radians=False)
        self.assertAlmostEqual(lon2, 45.0, places=8)
        self.assertAlmostEqual(lat2, 0.0, places=8)
        self.assertAlmostEqual(baz, 270.0, places=8)
        return

    def test_SphericalForward2(self):
        lon1 = 30.0
        lat1 = 0.0
        lon2, lat2, baz = crs.SphericalEarth.forward(lon1, lat1, 90.0, 5003778.767588614,
                                                     radians=False)
        self.assertAlmostEqual(lon2, 75.0, places=8)
        self.assertAlmostEqual(lat2, 0.0, places=8)
        self.assertAlmostEqual(baz, 270.0, places=8)
        return

    def test_SphericalForward3(self):
        lon1 = -120.0
        lat1 = 49.0
        lon2, lat2, baz = crs.SphericalEarth.forward(lon1, lat1, 310.0, 2000e3,
                                                     radians=False)
        self.assertAlmostEqual(lon2, -146.5118, places=4)
        self.assertAlmostEqual(lat2, 57.9998, places=4)
        self.assertAlmostEqual(baz, 108.48895148, places=4)
        return

    def test_SphericalInverse1(self):
        lon1 = 0.0
        lat1 = 0.0
        lon2 = -45.0
        lat2 = 0.0
        az, baz, dist = crs.SphericalEarth.inverse(lon1, lat1, lon2, lat2, radians=False)
        self.assertEqual(az, 270.0)
        self.assertEqual(baz, 90.0)
        self.assertAlmostEqual(dist, 5003778.767589, places=6)
        return

    def test_SphericalInverse2(self):
        lon1 = 32.0
        lat1 = -17.0
        lon2 = 38.0
        lat2 = 5.0
        az, baz, dist = crs.SphericalEarth.inverse(lon1, lat1, lon2, lat2, radians=False)
        self.assertAlmostEqual(az, 15.5977, places=4)
        self.assertAlmostEqual(baz, 194.9583, places=4)
        self.assertAlmostEqual(dist, 2533572.0748, places=2)
        return

    def test_SphericalInverse3(self):
        lon1 = 32.0
        lat1 = 5.0
        lon2 = 38.0
        lat2 = -17.0
        az, baz, dist = crs.SphericalEarth.inverse(lon1, lat1, lon2, lat2, radians=False)
        self.assertAlmostEqual(az, 165.0417, places=4)
        self.assertAlmostEqual(baz, 344.4023, places=4)
        self.assertAlmostEqual(dist, 2533572.0748, places=2)
        return

    def test_EllipsoidalEquatorialAzimuth(self):
        az, baz, _ = crs.LonLatWGS84.inverse(-40.0, 0.0, 55.0, 0.0)
        self.assertEqual(az, 90)
        self.assertEqual(baz, -90)

        az2, baz2, _ = crs.LonLatWGS84.inverse(180.0, 0.0, 5, 0.0)
        self.assertEqual(az2, -90)
        self.assertEqual(baz2, 90)
        return

    def test_EllipsoidalNearAntipodalInverse(self):
        az, baz, d = crs.LonLatWGS84.inverse(0.0, 30.0, 179.9, -29.9)
        az_, baz_, d_ = crs.LonLatWGS84_proj4.inverse(0.0, 30.0, 179.9, -29.9)
        self.assertAlmostEqual(az, az_, places=4)
        self.assertAlmostEqual(baz, baz_, places=4)
        self.assertAlmostEqual(d, d_, places=4)
        return

    def test_EllipsoidalForward(self):
        np.random.seed(43)
        for i in range(500):
            x = 360*np.random.rand() - 180
            y = 180*np.random.rand() - 90
            az = 360*np.random.rand() - 180
            d = 2e7*np.random.rand()
            x1, y1, baz = crs.LonLatWGS84.forward(x, y, az, d)
            x1_, y1_, baz_ = crs.LonLatWGS84_proj4.forward(x, y, az, d)

            self.assertAlmostEqual(x1, x1_, places=4)
            self.assertAlmostEqual(y1, y1_, places=4)
            self.assertAlmostEqual(baz, baz_, places=4)

    def test_EllipsoidalInverse(self):
        np.random.seed(43)
        for i in range(500):
            x1 = 360*np.random.rand() - 180
            y1 = 178*np.random.rand() - 89
            x2 = 360*np.random.rand() - 180
            y2 = 178*np.random.rand() - 89
            az, baz, d = crs.LonLatWGS84.inverse(x1, y1, x2, y2)
            az_, baz_, d_ = crs.LonLatWGS84_proj4.inverse(x1, y1, x2, y2)

            self.assertAlmostEqual(az, az_, places=4)
            self.assertAlmostEqual(baz, baz_, places=4)
            self.assertAlmostEqual(d, d_, places=2)

    def test_ConstructProj4(self):
        # Canonical constructor
        crs.Proj4CRS("+proj=longlat +datum=WGS84 +no_defs", "+ellps=WGS84")

        # Combine ellipsoid information into projection information
        crs.Proj4CRS("+proj=longlat +ellps=WGS84 +no_defs")

        # Fail if no geodetic information given
        self.assertRaises(karta.errors.CRSError,
                          lambda: crs.Proj4CRS("+proj=longlat +no_defs"))
        return

    def test_equal1(self):
        WGS84 = crs.Proj4CRS("+proj=longlat +datum=WGS84 +no_defs",
                              "+ellps=WGS84", name="WGS84 (Geographical)")
        WGS84_ = crs.Proj4CRS("+proj=longlat +datum=WGS84 +no_defs",
                               "+ellps=WGS84", name="WGS84 (Geographical)")
        self.assertTrue(WGS84 == WGS84_)
        return

    def test_equal2(self):
        WGS84 = crs.Proj4CRS("+proj=longlat +datum=WGS84 +no_defs",
                              "+ellps=WGS84", name="WGS84 (Geographical)")
        NAD83 = crs.Proj4CRS("+proj=longlat +datum=WGS84 +no_defs",
                              "+ellps=GRS80", name="NAD83 (Geographical)")
        self.assertTrue(not WGS84 == NAD83)
        return

    def test_not_equal1(self):
        WGS84 = crs.Proj4CRS("+proj=longlat +datum=WGS84 +no_defs",
                              "+ellps=WGS84", name="WGS84 (Geographical)")
        NAD83 = crs.Proj4CRS("+proj=longlat +datum=WGS84 +no_defs",
                              "+ellps=GRS80", name="NAD83 (Geographical)")
        self.assertTrue(WGS84 != NAD83)
        return

    def test_not_equal2(self):
        WGS84 = crs.Proj4CRS("+proj=longlat +datum=WGS84 +no_defs",
                              "+ellps=WGS84", name="WGS84 (Geographical)")

        WGS84_ = crs.Proj4CRS("+proj=longlat +datum=WGS84 +no_defs",
                               "+ellps=WGS84", name="WGS84 (Geographical)")
        self.assertTrue(not WGS84 != WGS84_)
        return

    def test_brent1(self):
        def forsythe(x):
            return x**3 - 2*x - 5
        self.assertAlmostEqual(2.094551482,
                geodesy.fzero_brent(2, 3, forsythe, 1e-12))
        return

    def test_brent2(self):
        self.assertAlmostEqual(0.7390851332,
                geodesy.fzero_brent(0, 1, lambda x: math.cos(x)-x, 1e-12))
        return

    def test_brent3(self):
        self.assertAlmostEqual(0.0,
                geodesy.fzero_brent(-1, 1, lambda x: math.sin(x)-x, 1e-12))
        return

    def test_brent4(self):
        self.assertAlmostEqual(0.0,
                geodesy.fzero_brent(0, 1, lambda x: math.sin(x)-x, 1e-12))
        return

    def test_brent_bracket_error(self):
        self.assertRaises(ValueError,
                geodesy.fzero_brent, 0.2, 2, lambda x: math.sin(x)-x, 1e-12)
        return


if __name__ == "__main__":
    unittest.main()

