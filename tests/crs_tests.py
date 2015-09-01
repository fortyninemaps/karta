
import unittest
import math
import karta.crs as crs

class TestCRS(unittest.TestCase):

    def test_get_proj4(self):
        c0 = crs.SphericalEarth
        self.assertEqual(c0.get_proj4(), "+proj=lonlat +a=6371009.000000 +b=6371009.000000 +no_defs")

        c1 = crs.LonLatWGS84
        self.assertEqual(c1.get_proj4(), "+proj=lonlat +a=6378137.0 +f=0.0033528106647474805")

        c2 = crs.NSIDCNorth
        self.assertEqual(c2.get_proj4(),
            "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 "
            "+units=m +datum=WGS84 +no_defs  +a=6378137.0 +f=0.0033528106647474805")
        return

    def test_get_wkt(self):
        c0 = crs.SphericalEarth
        self.assertEqual(c0.get_wkt(), '+proj=lonlat +a=6371009.000000 +b=6371009.000000 +no_defs')

        c1 = crs.LonLatWGS84
        self.assertEqual(c1.get_wkt(), 'GEOGCS["unnamed ellipse",DATUM["unknown",SPHEROID["unnamed",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]')

        c2 = crs.NSIDCNorth
        self.assertEqual(c2.get_wkt(),
            'PROJCS["unnamed",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9108"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",70],PARAMETER["central_meridian",0],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]')
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

    def test_ConstructProj4(self):
        # Canonical constructor
        crs.Proj4CRS("+proj=longlat +datum=WGS84 +no_defs", "+ellps=WGS84")

        # Combine ellipsoid information into projection information
        crs.Proj4CRS("+proj=longlat +ellps=WGS84 +no_defs")

        # Fail if no geodetic information given
        self.assertRaises(crs.CRSError,
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

if __name__ == "__main__":
    unittest.main()

