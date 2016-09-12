
import unittest
import math
import numpy as np
import karta.crs as crs
import karta.geodesy as geodesy
import karta.errors

class TestCRS(unittest.TestCase):

    def assertTuplesAlmostEqual(self, a, b, tol=1e-8):
        if len(a) != len(b):
            unittest.fail()
        for _a, _b in zip(a, b):
            if abs(_a-_b) > tol:
                unittest.fail()
        return

    def test_get_proj4(self):
        c0 = crs.SphericalEarth
        self.assertEqual(c0.get_proj4(), "+proj=lonlat +ellps=sphere +datum=WGS84")

        c1 = crs.LonLatWGS84
        proj4 = c1.get_proj4()
        self.assertTrue("+proj=lonlat" in proj4)
        self.assertTrue("+ellps=WGS84" in proj4)

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
        return

    def test_get_wkt(self):
        c0 = crs.SphericalEarth
        self.assertTrue(c0.wkt.startswith('GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]]'))

        c1 = crs.LonLatWGS84
        self.assertTrue('GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.25722' in c1.wkt)
        self.assertTrue('PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433' in c1.wkt)

        c2 = crs.NSIDCNorth
        self.assertTrue('PROJCS["unnamed",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.25722' in c2.wkt)
        #self.assertTrue('AUTHORITY["EPSG","7030"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329' in c2.wkt)
        #self.assertTrue('AUTHORITY["EPSG","9108"]],AUTHORITY["EPSG","4326"]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",70],PARAMETER["central_meridian",-45],PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]' in c2.wkt)
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
        d_ans = math.sqrt(a1**2 + a2**2)
        az_ans = geodesy.reduce_rad(1.5*math.pi - math.atan(a2/a1))
        baz_ans = geodesy.reduce_rad(az_ans - math.pi)
        az, baz, d = crs.Cartesian.inverse(lon0, lat0, lon1, lat1, radians=True)
        self.assertAlmostEqual(d_ans, d, places=12)
        self.assertAlmostEqual(az_ans, az, places=12)
        self.assertAlmostEqual(baz_ans, baz, places=12)
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
        self.assertAlmostEqual(lon2, 45.000084759104425, places=8)
        self.assertAlmostEqual(lat2, 0.0, places=8)
        self.assertAlmostEqual(baz, -90.0, places=8)
        return

    def test_SphericalForward2(self):
        lon1 = 30.0
        lat1 = 0.0
        lon2, lat2, baz = crs.SphericalEarth.forward(lon1, lat1, 90.0, 5003778.767588614,
                                                     radians=False)
        self.assertAlmostEqual(lon2, 75.00008475910442, places=8)
        self.assertAlmostEqual(lat2, 0.0, places=8)
        self.assertAlmostEqual(baz, -90.0, places=8)
        return

    def test_SphericalForward3(self):
        lon1 = -120.0
        lat1 = 49.0
        lon2, lat2, baz = crs.SphericalEarth.forward(lon1, lat1, 310.0, 2000e3,
                                                     radians=False)
        self.assertAlmostEqual(lon2, -146.51186194714958, places=6)
        self.assertAlmostEqual(lat2, 57.99979808258465, places=6)
        self.assertAlmostEqual(baz, 108.48890006687964, places=6)
        return

    def test_SphericalInverse1(self):
        lon1 = 0.0
        lat1 = 0.0
        lon2 = -45.0
        lat2 = 0.0
        az, baz, dist = crs.SphericalEarth.inverse(lon1, lat1, lon2, lat2, radians=False)
        self.assertEqual(az, -90.0)
        self.assertEqual(baz, 90.0)
        self.assertAlmostEqual(dist, 5003769.342810653, places=6)
        return

    def test_SphericalInverse2(self):
        lon1 = 32.0
        lat1 = -17.0
        lon2 = 38.0
        lat2 = 5.0
        az, baz, dist = crs.SphericalEarth.inverse(lon1, lat1, lon2, lat2, radians=False)
        self.assertAlmostEqual(az, 15.597740818516172, places=6)
        self.assertAlmostEqual(baz, -165.04174639642943, places=6)
        self.assertAlmostEqual(dist, 2533567.302751705, places=6)
        return

    def test_SphericalInverse3(self):
        lon1 = 32.0
        lat1 = 5.0
        lon2 = 38.0
        lat2 = -17.0
        az, baz, dist = crs.SphericalEarth.inverse(lon1, lat1, lon2, lat2, radians=False)
        self.assertAlmostEqual(az, 165.0417463964294, places=6)
        self.assertAlmostEqual(baz, -15.597740818516172, places=6)
        self.assertAlmostEqual(dist, 2533567.302751705, places=6)
        return

    def test_SphericalArea(self):
        r = 6378137.0
        x1 = 0.0
        x2 = 137.84490004377
        y1 = 40.0
        y2 = 41.79331020506
        S12 = karta.geodesy.spherical_area(r, x1, y1, x2, y2)
        self.assertAlmostEqual(abs(S12)/1e6, 84516702.1955, places=4)
        return

    def test_SphericalArea_dateline(self):
        r = 6378137.0
        x1 = 70.0
        x2 = 207.84490004377
        y1 = 40.0
        y2 = 41.79331020506
        S12 = karta.geodesy.spherical_area(r, x1, y1, x2, y2)
        self.assertAlmostEqual(abs(S12)/1e6, 84516702.1955, places=4)
        return


    def test_SphericalIntersection(self):
        ix = geodesy.intersection_spherical([(45, 10), (60, 10)],
                                            [(50, -10), (50, 20)])
        self.assertTuplesAlmostEqual(ix, (50, 10.075124337))

        ix = geodesy.intersection_spherical([(45, 0), (-140, 0)],
                                            [(50, -10), (50, 10)])
        self.assertTuplesAlmostEqual(ix, (50, 0))

        self.assertRaises(karta.errors.NoIntersection,
                lambda: geodesy.intersection_spherical([(45, 0), (230, 0)],
                                                       [(50, -10), (50, 10)]))
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
        az_, baz_, d_ = crs._LonLatWGS84.inverse(0.0, 30.0, 179.9, -29.9)
        self.assertAlmostEqual(az, az_, places=4)
        self.assertAlmostEqual(baz, baz_, places=4)
        self.assertAlmostEqual(d, d_, places=4)
        return

    def test_EllipsoidalForward(self):
        np.random.seed(43)
        for i in range(1000):
            x = 360*np.random.rand() - 180
            y = 180*np.random.rand() - 90
            az = 360*np.random.rand() - 180
            d = 2e7*np.random.rand()
            x1, y1, baz = crs.LonLatWGS84.forward(x, y, az, d)
            x1_, y1_, baz_ = crs._LonLatWGS84.forward(x, y, az, d)

            self.assertAlmostEqual(x1, x1_, places=4)
            self.assertAlmostEqual(y1, y1_, places=4)
            self.assertAlmostEqual(baz, baz_, places=4)

    def test_EllipsoidalInverse(self):
        np.random.seed(43)
        for i in range(1000):
            x1 = 360*np.random.rand() - 180
            y1 = 178*np.random.rand() - 89
            x2 = 360*np.random.rand() - 180
            y2 = 178*np.random.rand() - 89
            az, baz, d = crs.LonLatWGS84.inverse(x1, y1, x2, y2)
            az_, baz_, d_ = crs._LonLatWGS84.inverse(x1, y1, x2, y2)

            self.assertAlmostEqual(az, az_, places=4)
            self.assertAlmostEqual(baz, baz_, places=4)
            self.assertAlmostEqual(d, d_, places=2)

    def test_EllipsoidalArea(self):
        a = 6378137.0
        b = 6356752.314245
        x1 = 0.0
        x2 = 137.84490004377
        y1 = 40.0
        y2 = 41.79331020506
        S12 = karta.geodesy.ellipsoidal_area(a, b, x1, y1, x2, y2)
        self.assertAlmostEqual(abs(S12)/1e6, 84275623.42235, places=4)
        return

    def test_EllipsoidalArea_dateline(self):
        a = 6378137.0
        b = 6356752.314245
        x1 = 70.0
        x2 = 207.84490004377
        y1 = 40.0
        y2 = 41.79331020506
        S12 = karta.geodesy.ellipsoidal_area(a, b, x1, y1, x2, y2)
        self.assertAlmostEqual(abs(S12)/1e6, 84275623.42235, places=4)
        return

    def test_ConstructProj4(self):
        # Canonical constructor
        crs.ProjectedCRS("+proj=longlat +datum=WGS84 +no_defs", "+ellps=WGS84")

        # Combine ellipsoid information into projection information
        crs.ProjectedCRS("+proj=longlat +ellps=WGS84 +no_defs")

        # Fail if no geodetic information given
        self.assertRaises(karta.errors.CRSError,
                          lambda: crs.ProjectedCRS("+proj=longlat +no_defs"))
        return

    def test_equal1(self):
        WGS84 = crs.ProjectedCRS("+proj=longlat +datum=WGS84 +no_defs",
                                name="WGS84 (Geographical)")
        WGS84_ = crs.ProjectedCRS("+proj=longlat +datum=WGS84 +no_defs",
                                name="WGS84 (Geographical)")
        self.assertTrue(WGS84 == WGS84_)
        return

    def test_equal2(self):
        WGS84 = crs.ProjectedCRS("+proj=longlat +datum=WGS84 +no_defs",
                                name="WGS84 (Geographical)")
        NAD83 = crs.ProjectedCRS("+proj=longlat +datum=NAD83 +no_defs",
                                name="NAD83 (Geographical)")
        self.assertTrue(not WGS84 == NAD83)
        return

    def test_not_equal1(self):
        WGS84 = crs.ProjectedCRS("+proj=longlat +datum=WGS84 +no_defs",
                                name="WGS84 (Geographical)")
        NAD83 = crs.ProjectedCRS("+proj=longlat +datum=NAD83 +no_defs",
                                name="NAD83 (Geographical)")
        self.assertTrue(WGS84 != NAD83)
        return

    def test_not_equal2(self):
        WGS84 = crs.ProjectedCRS("+proj=longlat +datum=WGS84 +no_defs",
                                name="WGS84 (Geographical)")

        WGS84_ = crs.ProjectedCRS("+proj=longlat +datum=WGS84 +no_defs",
                                name="WGS84 (Geographical)")
        self.assertTrue(not WGS84 != WGS84_)
        return

    def test_export_wkt_same_as_parent(self):
        # Checks for bug that occured where shapefiles saved with NSIDCNorth
        # (WKT) were not loaded properly because some of the projection
        # information was mangled

        def proj4_asdict(s):
            keys = [a.split("=")[0] for a in s.split() if "=" in a]
            values = [a.split("=")[1] for a in s.split() if "=" in a]
            return dict(list(zip(keys, values)))

        wkt = crs.NSIDCNorth.get_wkt()
        new_crs = crs.crs_from_wkt(wkt)

        crsdict1 = proj4_asdict(crs.NSIDCNorth.get_proj4())
        crsdict2 = proj4_asdict(new_crs.get_proj4())
        for key in ("+proj", "+lat_0", "+lon_0", "+lat_ts", "+datum"):
            self.assertEqual(crsdict1[key], crsdict2[key])

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

    def test_datum_transform(self):
        lng, lat = crs.LonLatNAD27.transform(crs.LonLatNAD83, -107.5, 43.14)
        self.assertAlmostEqual(lng, -107.50062798611111, places=3)
        self.assertAlmostEqual(lat, 43.13996053333333, places=3)

class TestGeodesyFuncs(unittest.TestCase):

    def assertTuplesAlmostEqual(self, a, b, tol=1e-8):
        if len(a) != len(b):
            self.fail()
        for _a, _b in zip(a, b):
            if abs(_a-_b) > tol:
                self.fail()
        return

    def test_isbetween_circular(self):
        self.assertTrue(geodesy.isbetween_circular(90, 80, 100))
        self.assertTrue(geodesy.isbetween_circular(60, 90, -80))
        self.assertTrue(geodesy.isbetween_circular(50, 45, 224))
        self.assertFalse(geodesy.isbetween_circular(70, 90, 100))
        self.assertFalse(geodesy.isbetween_circular(60, 90, -110))
        self.assertFalse(geodesy.isbetween_circular(50, 45, 226))
        self.assertFalse(geodesy.isbetween_circular(50, 45, 225))
        return

    def test_cross_3d(self):
        self.assertEqual(geodesy.cross((0, 0, 1), (0, 0, 1)), (0, 0, 0))
        self.assertEqual(geodesy.cross((0, 0, 1), (1, 0, 0)), (0, 1, 0))
        self.assertEqual(geodesy.cross((1, 0, 0), (0, 0, 1)), (0, -1, 0))
        return

    def test_cart2sph(self):
        self.assertTuplesAlmostEqual(geodesy.cart2sph(1, 1, 1),
                                    (45.0, 35.2643896827))
        self.assertTuplesAlmostEqual(geodesy.cart2sph(1, 0, 1), (0.0, 45.0))
        self.assertTuplesAlmostEqual(geodesy.cart2sph(-1, 1, 0), (135.0, 0.0))
        return

if __name__ == "__main__":
    unittest.main()

