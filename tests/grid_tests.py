""" Unit tests for raster functions """

import unittest
import os
import numpy as np
import numpy.testing as npt
from test_helper import TESTDATA

import karta
from karta import RegularGrid

class RegularGridTests(unittest.TestCase):

    def setUp(self):
        pe = peaks(n=49)
        self.rast = RegularGrid((0.0, 0.0, 30.0, 30.0, 0.0, 0.0), values=pe)

    def test_get_resolution(self):
        grid = RegularGrid([0.0, 0.0, 25.0, 35.0, 10.0, 10.0])
        self.assertEqual(grid.resolution, (25.0, 35.0))

    def test_get_coordinates(self):
        cg = self.rast.coordinates()
        self.assertTrue(isinstance(cg, karta.raster.coordgen.CoordinateGenerator))
        self.assertEqual(cg.transform, self.rast.transform)

    def test_bbox(self):
        pe = peaks(n=49)
        pe[40:,:] = np.nan
        pe[:,:5] = np.nan
        grid = RegularGrid((0.0, 0.0, 30.0, 30.0, 0.0, 0.0), values=pe)
        self.assertEqual(grid.bbox, (0, 0, 30*49, 30*49))

    def test_data_bbox(self):
        pe = peaks(n=49)
        pe[40:,:] = np.nan
        pe[:,:5] = np.nan
        grid = RegularGrid((0.0, 0.0, 30.0, 30.0, 0.0, 0.0), values=pe)
        self.assertEqual(grid.data_bbox, (5*30, 0, 30*49, 30*40))

    def test_add_rgrid(self):
        rast2 = RegularGrid(self.rast.transform,
                            values=np.random.random(self.rast.size))
        res = self.rast + rast2
        self.assertTrue(np.all(res[:,:] == self.rast[:,:]+rast2[:,:]))

    def test_sub_rgrid(self):
        rast2 = RegularGrid(self.rast.transform,
                            values=np.random.random(self.rast.size))
        res = self.rast - rast2
        self.assertTrue(np.all(res[:,:] == self.rast[:,:]-rast2[:,:]))

    def test_center_coords(self):
        grid = RegularGrid((0.0, 0.0, 30.0, 30.0, 0.0, 0.0),
                           values=np.zeros([49, 49]))
        ans = np.meshgrid(np.arange(15.0, 1471.0, 30.0),
                          np.arange(15.0, 1471.0, 30.0))
        self.assertEqual(0.0, np.sum(grid.center_coords()[0] - ans[0]))
        self.assertEqual(0.0, np.sum(grid.center_coords()[1] - ans[1]))

    def test_center_coords_skewed(self):
        grid = RegularGrid((15.0, 15.0, 30.0, 30.0, 20.0, 10.0),
                           values=np.zeros([5, 5]))
        X, Y = grid.center_coords()
        self.assertEqual(X[0,0], 40.0)
        self.assertEqual(Y[0,0], 35.0)
        self.assertEqual(X[-1,0], 120.0)
        self.assertEqual(Y[-1,0], 155.0)
        self.assertEqual(X[-1,-1], 240.0)
        self.assertEqual(Y[-1,-1], 195.0)

    def test_aschunks(self):
        grid = RegularGrid((0, 0, 1, 1, 0, 0), values=peaks(512))
        for chunk in grid.aschunks(size=(64, 64)):
            self.assertEqual(chunk.size, (64, 64))

    def test_aschunks_overlap(self):
        grid = RegularGrid((0, 0, 1, 1, 0, 0), values=peaks(512))
        count = 0
        for chunk in grid.aschunks(size=(64, 64), overlap=(32, 32)):
            count += 1
        self.assertEqual(count, 256)

    def test_index_multiband(self):
        p = peaks(512)
        grid = RegularGrid((0, 0, 1, 1, 0, 0),
                           values=np.dstack([p,
                                             np.where(p > 0, p, np.nan),
                                             p]),
                            nodata_value = np.nan)
        grid[grid.data_mask_full]

    def test_apply(self):
        msk = np.zeros([8, 8], dtype=np.bool)
        msk[3, 2] = True
        msk[:2,3:] = True
        val = np.arange(64, dtype=np.float64).reshape([8,8])
        val_ = val.copy()
        val[msk] = -1
        grid = RegularGrid([0, 0, 1, 1, 0, 0], values=val, nodata_value=-1)
        newgrid = grid.apply(lambda x: x**2)

        npt.assert_equal(newgrid[msk], -1)
        npt.assert_equal(newgrid[:,:,0][~msk], val_[~msk]**2)

    def test_apply_inplace(self):
        msk = np.zeros([8, 8], dtype=np.bool)
        msk[3, 2] = True
        msk[:2,3:] = True
        val = np.arange(64, dtype=np.float64).reshape([8,8])
        val_ = val.copy()
        val[msk] = -1
        grid = RegularGrid([0, 0, 1, 1, 0, 0], values=val, nodata_value=-1)
        grid.apply(lambda x: x**2, inplace=True)

        self.assertTrue(np.all(grid[:,:,0][msk] == -1))
        self.assertTrue(np.all(grid[:,:,0][~msk] == val_[~msk]**2))

    def test_merge(self):
        grid1 = RegularGrid([10, 20, 1, 1, 0, 0], values=np.ones([8, 8]))
        grid2 = RegularGrid([7, 22, 1, 1, 0, 0], values=2*np.ones([4, 6]))
        grid3 = RegularGrid([12, 15, 1, 1, 0, 0], values=3*np.ones([5, 5]))
        grid_combined = karta.raster.merge([grid1, grid2, grid3])
        self.assertEqual(grid_combined.transform, (7.0, 15.0, 1.0, 1.0, 0.0, 0.0))
        self.assertEqual(grid_combined.size, (13, 11))
        self.assertEqual(np.sum(np.isnan(grid_combined[:,:])), 42)

    def test_merge_weighted(self):
        grid1 = RegularGrid([10, 20, 1, 1, 0, 0], values=np.ones([8, 8]))
        grid2 = RegularGrid([7, 22, 1, 1, 0, 0], values=2*np.ones([4, 6]))
        grid3 = RegularGrid([12, 19, 1, 1, 0, 0], values=3*np.ones([5, 5]))
        grid_combined = karta.raster.merge([grid1, grid2, grid3], weights=[1, 2, 3])
        self.assertAlmostEqual(grid_combined[4,4,0], 1.66666666666)
        self.assertAlmostEqual(grid_combined[2,8,0], 2.5)
        self.assertAlmostEqual(grid_combined[4,5,0], 2.33333333333)

    def test_merge_multiband(self):
        grid3a = RegularGrid([0, 0, 1, 1, 0, 0],
                             values=np.array([1,2,3]) * np.ones((16, 16, 3)))
        grid3b = RegularGrid([4, 4, 1, 1, 0, 0],
                             values=np.array([2,3,4]) * np.ones((16, 16, 3)))
        grid3_mosaic = karta.raster.merge([grid3a, grid3b])
        self.assertEqual(np.nansum(grid3_mosaic[:,:,0]), 552)
        self.assertEqual(np.nansum(grid3_mosaic[:,:,1]), 920)
        self.assertEqual(np.nansum(grid3_mosaic[:,:,2]), 1288)

    def test_align_origin(self):
        xx, yy = np.meshgrid(np.linspace(0, 1, 50), np.linspace(0, 1, 30))
        zz = 2.0*xx**2 - 3.0*yy**2
        grid = RegularGrid((27, 53, 5, 5, 0, 0), values=zz)
        new_grid = grid._align_origin(5, 5, method='linear')
        self.assertEqual(new_grid.origin, (25, 55))
        self.assertTrue(np.isnan(new_grid.values[0, 0]))

    def test_resample_nearest(self):
        # use linear function so that nearest neighbour and linear interp are
        # exact
        def makegrid(start, finish, n, res):
            xx, yy = np.meshgrid(np.linspace(start, finish, n),
                                 np.linspace(start, finish, n))
            zz = 2.0*xx - 3.0*yy
            return RegularGrid((0.0, 0.0, res, res, 0.0, 0.0), values=zz)

        # node numbers from a line with extreme edges at [0, 1]
        g = makegrid(0.0, 1.0-2.0/300, 150, 2.0)
        sol = makegrid(0.0, 1.0-6.0/300, 50, 6.0)
        gnew = g.resample(6.0, 6.0, method='nearest')
        residue = gnew[:,:] - sol[:,:]
        self.assertTrue(np.max(np.abs(residue)) < 1e-12)

    def test_resample_linear(self):
        # use linear function so that nearest neighbour and linear interp are
        # exact
        def makegrid(start, finish, n, res):
            xx, yy = np.meshgrid(np.linspace(start, finish, n),
                                 np.linspace(start, finish, n))
            zz = 2.0*xx - 3.0*yy
            return RegularGrid((0.0, 0.0, res, res, 0.0, 0.0), values=zz)

        # node numbers from a line with extreme edges at [0, 1]
        g = makegrid(0.0, 1.0-2.0/300, 150, 2.0)
        sol = makegrid(0.0, 1.0-6.0/300, 50, 6.0)
        gnew = g.resample(6.0, 6.0, method='linear')
        residue = gnew[:,:] - sol[:,:]
        self.assertTrue(np.max(np.abs(residue)) < 1e-12)

    def test_sample_nearest_out_of_bounds(self):
        g = RegularGrid([0, 0, 1, 1, 0, 0], values=np.ones((10, 10)))
        v = g.sample_nearest(np.array([7, 9, 12, 15]), np.array([3, 1, -1, 1]))
        self.assertEqual(v[0][0], 1.0)
        self.assertEqual(v[0][1], 1.0)
        self.assertTrue(np.isnan(v[0][2]))
        self.assertTrue(np.isnan(v[0][3]))

    def test_sample_linear_out_of_bounds(self):
        g = RegularGrid([0, 0, 1, 1, 0, 0], values=np.ones((10, 10)))
        v = g.sample_bilinear(np.array([7, 9, 12, 15]), np.array([3, 1, -1, 1]))
        self.assertEqual(v[0][0], 1.0)
        self.assertEqual(v[0][1], 1.0)
        self.assertTrue(np.isnan(v[0][2]))
        self.assertTrue(np.isnan(v[0][3]))

    def test_resample_multiband(self):
        grid = RegularGrid((0, 0, 1, 1, 0, 0),
                           values=np.dstack([np.ones((64, 64)),
                                             2*np.ones((64, 64)),
                                             3*np.ones((64, 64))]))
        grid2 = grid.resample(0.5, 0.5)
        self.assertEqual(grid2[0,0,0], 1.0)
        self.assertEqual(grid2[0,0,1], 2.0)
        self.assertEqual(grid2[0,0,2], 3.0)

    def test_sample_nearest(self):
        grid = RegularGrid([0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
                           values=np.array([[0, 1], [1, 0.5]]))
        self.assertEqual(grid.sample_nearest(0.6, 0.7), 0.0)
        self.assertEqual(grid.sample_nearest(0.6, 1.3), 1.0)
        self.assertEqual(grid.sample_nearest(1.4, 0.3), 1.0)
        self.assertEqual(grid.sample_nearest(1.6, 1.3), 0.5)

    def test_sample_nearest_vector(self):
        grid = RegularGrid([0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
                           values=np.arange(64).reshape([8,8]))
        res = grid.sample_nearest(np.arange(1,7,0.5), np.arange(2,5,0.25))
        self.assertEqual(res.shape, (1, 12))

    def test_sample_nearest_array(self):
        grid = RegularGrid([0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
                           values=np.arange(64).reshape([8,8]))
        X, Y = np.meshgrid(np.arange(1, 7, 0.5), np.arange(2, 5, 0.25))
        res = grid.sample_nearest(X, Y)
        self.assertEqual(res.shape, (1, 12, 12))

    def test_sample_nearest_array_order(self):
        grid = RegularGrid((0, 0, 1, 1, 0, 0),
                           values=np.dstack([np.ones((64, 64)),
                                             2*np.ones((64,64)),
                                             3*np.ones((64,64))]))

        X, Y = np.meshgrid(np.arange(2, 10), np.arange(5, 13))
        val = grid.sample_nearest(X, Y)
        self.assertEqual(val[0,0,0], 1.0)
        self.assertEqual(val[1,0,0], 2.0)
        self.assertEqual(val[2,0,0], 3.0)

    def test_sample_nearest_skewed(self):
        grid = RegularGrid([0.0, 0.0, 1.0, 1.0, 0.5, 0.2],
                           values=np.array([[0, 1], [1, 0.5]]))
        self.assertEqual(grid.sample_nearest(1, 0.75), 0.0)
        self.assertEqual(grid.sample_nearest(1.5, 1.05), 1.0)
        self.assertEqual(grid.sample_nearest(1.2, 1.4), 1.0)
        self.assertEqual(grid.sample_nearest(2.0, 1.7), 0.5)

    def test_sample_bilinear(self):
        grid = RegularGrid([0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
                           values=np.array([[0, 1], [1, 0.5]]))
        self.assertEqual(grid.sample_bilinear(1.0, 1.0), 0.625)

    def test_sample_bilinear_vector(self):
        grid = RegularGrid([0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
                           values=np.arange(64, dtype=np.float64).reshape([8,8]))
        res = grid.sample_bilinear(np.arange(1,7,0.5), np.arange(2,5,0.25))
        self.assertEqual(res.shape, (1, 12))

    def test_sample_bilinear_array(self):
        grid = RegularGrid([0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
                           values=np.arange(64, dtype=np.float64).reshape([8,8]))
        X, Y = np.meshgrid(np.arange(1, 7, 0.5), np.arange(2, 5, 0.25))
        res = grid.sample_bilinear(X, Y)
        self.assertEqual(res.shape, (1, 12, 12))


    def test_sample_bilinear_int(self):
        grid = RegularGrid([0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
                           values=np.array([[0, 2], [2, 1]], dtype=np.int32))
        self.assertEqual(grid.sample_bilinear(1.0, 1.0), 1)

    def test_sample_bilinear_uint(self):
        grid = RegularGrid([0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
                           values=np.array([[0, 2], [2, 1]], dtype=np.uint16))
        self.assertEqual(grid.sample_bilinear(1.0, 1.0), 1)

    def test_sample_bilinear_uint8(self):
        grid = RegularGrid([0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
                           values=np.array([[0, 2], [2, 1]], dtype=np.uint8))
        self.assertEqual(grid.sample_bilinear(1.0, 1.0), 1)

    def test_sample_bilinear_multiband(self):
        grid = RegularGrid([0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
                           values=np.dstack([[[0, 1], [1, 0.5]], [[1, 2], [2, 1.5]]]))
        res = grid.sample_bilinear(1.0, 1.0)
        self.assertTrue(np.all(res == np.array([0.625, 1.625])))

    def test_sample_bilinear_skewed(self):
        grid = RegularGrid([0.0, 0.0, 1.0, 1.0, 0.5, 0.2],
                           values=np.array([[0, 1], [1, 0.5]]))
        self.assertEqual(grid.sample_bilinear(1.5, 1.2), 0.625)

    def test_sample_bilinear2(self):
        grid = RegularGrid([0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
                           values=np.array([[0, 1], [1, 0.5]]))
        xi, yi = np.meshgrid(np.linspace(0.5, 1.5), np.linspace(0.5, 1.5))
        z = grid.sample_bilinear(xi, yi).ravel()
        self.assertEqual([z[400], z[1200], z[1550], z[2120]],
                         [0.16326530612244894, 0.48979591836734693,
                          0.63265306122448983, 0.74052478134110788])

    def test_sample_multipoint(self):
        grid = RegularGrid([0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
                           values=np.array([[0, 1], [1, 0.5]]))
        mp = karta.Multipoint([(0.6, 0.7), (0.6, 1.3), (1.4, 0.3), (1.6, 1.3)],
                              crs=grid.crs)
        self.assertTrue(np.all(np.allclose(grid.sample(mp, method="nearest"),
                                           np.array([0.0, 1.0, 1.0, 0.5]))))

    def test_sample_2d_array(self):
        grid = RegularGrid([0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
                           values=np.array([[0, 1], [1, 0.5]]))
        xi = np.array([[0.7, 1.3], [0.7, 1.3]])
        yi = np.array([[0.7, 0.7], [1.3, 1.3]])
        z = grid.sample(xi, yi, method="nearest")
        self.assertTrue(np.all(z == np.array([[0, 1], [1, 0.5]])))

    def test_vertex_coords(self):
        grid = RegularGrid((0.0, 0.0, 30.0, 30.0, 0.0, 0.0), values=np.zeros([49, 49]))
        ans = np.meshgrid(np.arange(15.0, 1486.0, 30.0),
                          np.arange(15.0, 1486.0, 30.0))
        self.assertTrue(np.sum(grid.vertex_coords()[0] - ans[0]) < 1e-10)
        self.assertTrue(np.sum(grid.vertex_coords()[1] - ans[1]) < 1e-10)

    def test_vertex_coords_skewed(self):
        grid = RegularGrid((0.0, 0.0, 30.0, 30.0, 20.0, 10.0), values=np.zeros([5, 5]))
        ans = np.meshgrid(np.arange(15.0, 1486.0, 30.0),
                          np.arange(15.0, 1486.0, 30.0))
        self.assertTrue(np.sum(self.rast.vertex_coords()[0] - ans[0]) < 1e-10)
        self.assertTrue(np.sum(self.rast.vertex_coords()[1] - ans[1]) < 1e-10)

    def test_get_extent(self):
        ereg = (0.0, 1470.0, 0.0, 1470.0)
        creg = (15.0, 1455.0, 15.0, 1455.0)
        self.assertEqual(self.rast.get_extent(reference='center'), creg)
        self.assertEqual(self.rast.get_extent(reference='edge'), ereg)

    def test_get_extent_crs(self):
        pe = peaks(n=49)
        crs = karta.crs.ProjectedCRS("+proj=utm +zone=12 +ellps=WGS84 +north=True", "UTM 12N (WGS 84)")
        rast_utm12N = RegularGrid((0.0, 0.0, 10000.0, 10000.0, 0.0, 0.0),
                                        values=pe,
                                        crs=crs)
        a,b,c,d = rast_utm12N.get_extent(reference='center',
                                         crs=karta.crs.LonLatWGS84)
        self.assertAlmostEqual(a, -115.45687156)
        self.assertAlmostEqual(b, -111.13480112)
        self.assertAlmostEqual(c, 0.0450996517)
        self.assertAlmostEqual(d, 4.3878488543)

    def test_get_data_extent(self):
        grid = RegularGrid((0, 0, 1, 1, 0, 0), values=np.zeros((128, 128), dtype=np.float64))
        grid[:,:4] = grid.nodata
        grid[:,-7:] = grid.nodata
        grid[:21,:] = grid.nodata
        grid[-17:,:] = grid.nodata
        self.assertEqual(grid.get_data_extent("edge"), (4, 121, 21, 111))
        self.assertEqual(grid.get_data_extent("center"), (4.5, 120.5, 21.5, 110.5))

    def test_minmax_nodata(self):
        values = np.array([[4, 5, 3], [4, 2, -9], [3, 6, 1]])

        self.rast = RegularGrid((0.0, 0.0, 30.0, 30.0, 0.0, 0.0),
                                values=values, nodata_value=-9)
        minmax = self.rast.minmax()
        self.assertEqual(minmax, (1, 6))

    def test_minmax_nodata2(self):
        values = -9*np.ones([3,3])

        self.rast = RegularGrid((0.0, 0.0, 30.0, 30.0, 0.0, 0.0),
                                values=values, nodata_value=-9)
        minmax = self.rast.minmax()
        self.assertTrue(np.isnan(minmax[0]))
        self.assertTrue(np.isnan(minmax[1]))

    def test_minmax(self):
        mx = self.rast.max()
        self.assertEqual(mx, 8.075173545159231)

        mn = self.rast.min()
        self.assertEqual(mn, -6.5466445243204294)

        minmax = self.rast.minmax()
        self.assertEqual(minmax, (-6.5466445243204294, 8.075173545159231))

    def test_clip(self):
        clipped = self.rast.clip(500, 950, 500, 950)
        self.assertEqual(clipped.size, (15, 15))
        self.assertEqual(clipped.transform, (510, 510, 30, 30, 0, 0))
        X, Y = clipped.center_coords()
        self.assertEqual(X[0,0], 525)
        self.assertEqual(X[0,-1], 945)
        self.assertEqual(Y[0,0], 525)
        self.assertEqual(Y[-1,0], 945)

    def test_clip_to_extent(self):
        proto = RegularGrid((500, 500, 30, 30, 0, 0), np.zeros((15,15)))
        clipped = self.rast.clip(*proto.get_extent("edge"))
        self.assertEqual(clipped.size, (15, 15))
        self.assertEqual(clipped.transform, (510, 510, 30, 30, 0, 0))
        X, Y = clipped.center_coords()
        self.assertEqual(X[0,0], 525)
        self.assertEqual(X[0,-1], 945)
        self.assertEqual(Y[0,0], 525)
        self.assertEqual(Y[-1,0], 945)

    def test_resize_smaller(self):
        proto = RegularGrid((500, 500, 30, 30, 0, 0), values=peaks(50))
        newgrid = proto.resize([620, 650, 1370, 1310])
        self.assertEqual(newgrid.transform, (620.0, 650.0, 30.0, 30.0, 0.0, 0.0))
        self.assertTrue(np.all(newgrid[:,:] == proto[5:27,4:29]))

    def test_resize_larger(self):
        proto = RegularGrid((500, 500, 30, 30, 0, 0), values=peaks(50))
        newgrid = proto.resize([380, 320, 380+30*60, 320+30*62])
        self.assertEqual(newgrid.transform, (380.0, 320.0, 30.0, 30.0, 0.0, 0.0))
        self.assertTrue(np.all(newgrid[6:56,4:54] == proto[:,:]))
        self.assertTrue(np.isnan(newgrid[0,0]))

    def test_resize_lower_left(self):
        proto = RegularGrid((500, 500, 30, 30, 0, 0), values=peaks(50))
        newgrid = proto.resize([380, 320, 380+30*30, 320+30*32])
        self.assertEqual(newgrid.transform, (380.0, 320.0, 30.0, 30.0, 0.0, 0.0))
        self.assertTrue(np.all(newgrid[6:,4:] == proto[:26,:26]))

    def test_resize_upper_right(self):
        proto = RegularGrid((500, 500, 30, 30, 0, 0), values=peaks(50))
        newgrid = proto.resize([1940, 1910, 1940+30*10, 1910+30*7])
        self.assertEqual(newgrid.transform, (1940.0, 1910.0, 30.0, 30.0, 0.0, 0.0))
        self.assertTrue(np.all(newgrid[:3,:2] == proto[-3:,-2:]))

    def test_data_mask_nan(self):
        T = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0]
        v = np.arange(64, dtype=np.float64).reshape([8, 8])
        v[0,2:7] = np.nan
        g = RegularGrid(T, values=v, nodata_value=np.nan)
        self.assertEqual(np.sum(g.data_mask), 59)

    def test_data_mask_nonnan(self):
        T = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0]
        v = np.arange(64, dtype=np.int8).reshape([8, 8])
        v[0,2:7] = -1
        g = RegularGrid(T, values=v, nodata_value=-1)
        self.assertEqual(np.sum(g.data_mask), 59)

    def test_mask_poly(self):
        t = -np.linspace(0, 2*np.pi, 200)
        xp = ((2+np.cos(7*t)) * np.cos(t+0.3) + 4) * 12
        yp = ((2+np.cos(7*t)) * np.sin(t+0.2) + 4) * 12
        poly = karta.Polygon(zip(xp, yp), crs=karta.crs.Cartesian)
        grid = RegularGrid([0.0, 0.0, 0.1, 0.1, 0.0, 0.0],
                           values=np.arange(1e6).reshape(1000, 1000),
                           crs=karta.crs.Cartesian)
        masked_grid = grid.mask_by_poly(poly)
        self.assertEqual(int(np.nansum(masked_grid[:,:])), 97048730546)

    def test_mask_poly_inplace(self):
        t = -np.linspace(0, 2*np.pi, 200)
        xp = ((2+np.cos(7*t)) * np.cos(t+0.3) + 4) * 12
        yp = ((2+np.cos(7*t)) * np.sin(t+0.2) + 4) * 12
        poly = karta.Polygon(zip(xp, yp), crs=karta.crs.Cartesian)
        grid = RegularGrid([0.0, 0.0, 0.1, 0.1, 0.0, 0.0],
                           values=np.arange(1e6).reshape(1000, 1000),
                           crs=karta.crs.Cartesian)
        grid.mask_by_poly(poly, inplace=True)
        self.assertEqual(int(np.nansum(grid[:,:])), 97048730546)

    def test_mask_poly_partial(self):
        t = -np.linspace(0, 2*np.pi, 200)
        xp = ((2+np.cos(7*t)) * np.cos(t+0.3) + 2) * 12
        yp = ((2+np.cos(7*t)) * np.sin(t+0.2) + 2) * 12
        poly = karta.Polygon(zip(xp, yp), crs=karta.crs.Cartesian)
        grid = RegularGrid([0.0, 0.0, 0.1, 0.1, 0.0, 0.0],
                           values=np.arange(1e6).reshape(1000, 1000),
                           crs=karta.crs.Cartesian)
        masked_grid = grid.mask_by_poly(poly)
        self.assertEqual(masked_grid.data_mask.sum(), 181424)

    def test_mask_poly_partial2(self):
        # this case has the first transect begin off-grid, which has caused
        # problems in the past
        g = RegularGrid([0, 0, 1, 1, 0, 0], values=np.ones((7, 7)))
        p = karta.Polygon([(-2, 3), (8, -5), (8, -1), (-2, 7)])
        gc = g.mask_by_poly(p)
        self.assertEqual(gc.data_mask.sum(), 20)

    def test_mask_poly_multiple(self):
        # mask by multiple polygons
        t = -np.linspace(0, 2*np.pi, 200)
        xp = ((2+np.cos(7*t)) * np.cos(t+0.3) + 4) * 4 + 15
        yp = ((2+np.cos(7*t)) * np.sin(t+0.2) + 4) * 4 + 72
        poly = karta.Polygon(zip(xp, yp), crs=karta.crs.Cartesian)
        xp2 = ((2+np.cos(7*t)) * np.cos(t+0.3) + 4) * 6 + 40
        yp2 = ((2+np.cos(7*t)) * np.sin(t+0.2) + 4) * 6 + 30
        poly2 = karta.Polygon(zip(xp2, yp2), crs=karta.crs.Cartesian)
        grid = RegularGrid([0.0, 0.0, 0.1, 0.1, 0.0, 0.0],
                           values=np.arange(1e6).reshape(1000, 1000),
                           crs=karta.crs.Cartesian)
        masked_grid = grid.mask_by_poly([poly, poly2])
        self.assertEqual(int(np.nansum(masked_grid[:,:])), 47081206720)

    def test_mask_poly_multiband(self):
        t = -np.linspace(0, 2*np.pi, 200)
        xp = ((2+np.cos(7*t)) * np.cos(t+0.3) + 4) * 12
        yp = ((2+np.cos(7*t)) * np.sin(t+0.2) + 4) * 12
        poly = karta.Polygon(zip(xp, yp), crs=karta.crs.Cartesian)
        grid = RegularGrid([0.0, 0.0, 0.1, 0.1, 0.0, 0.0],
                           values=np.broadcast_to(
                               np.atleast_3d(np.arange(1e6).reshape(1000, 1000)),
                               (1000, 1000, 3)),
                           crs=karta.crs.Cartesian)
        masked_grid = grid.mask_by_poly(poly)
        self.assertEqual(int(np.nansum(masked_grid[:,:])), 97048730546*3)

    def test_mask_poly_inverted(self):
        # grids where transform dy is negative (common for geotiffs)
        t = -np.linspace(0, 2*np.pi, 200)
        xp = ((2+np.cos(7*t)) * np.cos(t+0.3) + 4) * 12
        yp = ((2+np.cos(7*t)) * np.sin(t+0.2) + 4) * 12
        poly = karta.Polygon(zip(xp, yp), crs=karta.crs.Cartesian)
        grid = RegularGrid([0.0, 100, 0.1, -0.1, 0.0, 0.0],
                           values=np.arange(1e6).reshape(1000, 1000),
                           crs=karta.crs.Cartesian)
        masked_grid = grid.mask_by_poly(poly)
        self.assertEqual(int(np.nansum(masked_grid[:,:])), 97048730546)

    def test_mask_multipoly(self):
        t = -np.linspace(0, 2*np.pi, 200)

        coords = []
        for dx, dy in [(60, 30), (45, 80), (25, 35)]:
            xp = (2+np.cos(7*t)) * np.cos(t+0.3) * 6 + dx
            yp = (2+np.cos(7*t)) * np.sin(t+0.2) * 6 + dy
            coords.append([list(zip(xp, yp))])

        poly = karta.Multipolygon(coords, crs=karta.crs.Cartesian)
        grid = RegularGrid([0.0, 0, 0.1, 0.1, 0.0, 0.0],
                           values=np.arange(1e6).reshape(1000, 1000),
                           crs=karta.crs.Cartesian)
        masked_grid = grid.mask_by_poly(poly)

        self.assertEqual(int(np.nansum(masked_grid[:,:])), 73399874364)

    def test_get_positions(self):
        grid = RegularGrid([0.0, 0.0, 1.0, 1.0, 0.0, 0.0],
                                 values=np.zeros((3,3)))
        (i, j) = grid.get_positions(1.5, 1.5)
        self.assertEqual((i,j), (1.0, 1.0))
        (i, j) = grid.get_positions(2.0, 1.5)
        self.assertEqual((i,j), (1.0, 1.5))

    def test_get_indices(self):
        ind = self.rast.get_indices(15.0, 15.0)
        self.assertEqual(tuple(ind), (0, 0))
        ind = self.rast.get_indices(1455.0, 1455.0)
        self.assertEqual(tuple(ind), (48, 48))

    def test_get_indices_onevec(self):
        ind = self.rast.get_indices([15.0], [15.0])
        self.assertEqual(tuple(ind), (0, 0))
        ind = self.rast.get_indices(1455.0, 1455.0)
        self.assertEqual(tuple(ind), (48, 48))

    def test_get_indices_vec(self):
        ind = self.rast.get_indices(np.arange(15.0, 1470, 5),
                                    np.arange(15.0, 1470, 5))

        xi = np.array([ 0,  0,  0,  0,  1,  1,  1,  1,  1,  2,  2,  2,  2, 2,
            2, 2,  3, 3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  4,  5,  5, 5,  5,
            5,  6, 6,  6,  6,  6,  6,  6,  7,  7,  7,  7,  7,  8,  8, 8,  8, 8,
            8, 8,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11,
            11, 11, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 14, 14, 14,
            14, 14, 14, 14, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 17,
            17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 20,
            20, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 22, 22, 22, 22, 22, 22,
            22, 23, 23, 23, 23, 23, 24, 24, 24, 24, 24, 24, 24, 25, 25, 25, 25,
            25, 26, 26, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, 28, 28, 28, 28,
            28, 28, 28, 29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 30, 31, 31,
            31, 31, 31, 32, 32, 32, 32, 32, 32, 32, 33, 33, 33, 33, 33, 34, 34,
            34, 34, 34, 34, 34, 35, 35, 35, 35, 35, 36, 36, 36, 36, 36, 36, 36,
            37, 37, 37, 37, 37, 38, 38, 38, 38, 38, 38, 38, 39, 39, 39, 39, 39,
            40, 40, 40, 40, 40, 40, 40, 41, 41, 41, 41, 41, 42, 42, 42, 42, 42,
            42, 42, 43, 43, 43, 43, 43, 44, 44, 44, 44, 44, 44, 44, 45, 45, 45,
            45, 45, 46, 46, 46, 46, 46, 46, 46, 47, 47, 47, 47, 47, 48, 48, 48,
            48, 48, 48])


        yi = np.array([ 0,  0,  0,  0,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,
            2, 2,  3, 3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  4,  5,  5,  5,
            5, 5,  6, 6,  6,  6,  6,  6,  6,  7,  7,  7,  7,  7,  8,  8,  8,
            8, 8,  8, 8,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10, 10, 10, 11,
            11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 14,
            14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16,
            16, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19,
            19, 20, 20, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 22, 22, 22, 22,
            22, 22, 22, 23, 23, 23, 23, 23, 24, 24, 24, 24, 24, 24, 24, 25, 25,
            25, 25, 25, 26, 26, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, 28, 28,
            28, 28, 28, 28, 28, 29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 30,
            31, 31, 31, 31, 31, 32, 32, 32, 32, 32, 32, 32, 33, 33, 33, 33, 33,
            34, 34, 34, 34, 34, 34, 34, 35, 35, 35, 35, 35, 36, 36, 36, 36, 36,
            36, 36, 37, 37, 37, 37, 37, 38, 38, 38, 38, 38, 38, 38, 39, 39, 39,
            39, 39, 40, 40, 40, 40, 40, 40, 40, 41, 41, 41, 41, 41, 42, 42, 42,
            42, 42, 42, 42, 43, 43, 43, 43, 43, 44, 44, 44, 44, 44, 44, 44, 45,
            45, 45, 45, 45, 46, 46, 46, 46, 46, 46, 46, 47, 47, 47, 47, 47, 48,
            48, 48, 48, 48, 48])

        self.assertTrue(np.all(ind[0] == xi))
        self.assertTrue(np.all(ind[1] == yi))

    def test_profile(self):
        path = karta.Line([(15.0, 15.0), (1484.0, 1484.0)], crs=karta.crs.Cartesian)
        pts, z = self.rast.profile(path, resolution=42.426406871192853, method="nearest")
        expected = self.rast[:,:].diagonal()
        self.assertEqual(len(pts), 49)
        self.assertTrue(np.allclose(z, expected))

    def test_gridpoints_float64(self):
        # should use fast C version
        np.random.seed(49)
        x = np.random.rand(20000)*10.0-5.0
        y = np.random.rand(20000)*10.0-5.0
        z = x**2+y**3

        T = [-5.0, -5.0, 0.25, 0.25, 0.0, 0.0]
        grid = karta.raster.gridpoints(x, y, z, T, karta.crs.Cartesian)

        Xg, Yg = grid.center_coords()
        self.assertTrue(np.sum(np.abs(Xg**2+Yg**3-grid[:,:,0]))/Xg.size < 0.45)

    def test_gridpoints_float32(self):
        # should use Python fallback
        np.random.seed(49)
        x = np.random.rand(20000).astype(np.float32)*10.0-5.0
        y = np.random.rand(20000).astype(np.float32)*10.0-5.0
        z = x**2+y**3

        T = [-5.0, -5.0, 0.25, 0.25, 0.0, 0.0]
        grid = karta.raster.gridpoints(x, y, z, T, karta.crs.Cartesian)

        Xg, Yg = grid.center_coords()
        self.assertTrue(np.sum(np.abs(Xg**2+Yg**3-grid[:,:,0]))/Xg.size < 0.45)

    def test_set_nodata(self):
        v = np.arange(64, dtype=np.float64).reshape([8,8])
        v[2:4, 5:7] = -1
        grid = RegularGrid([0, 0, 1, 1, 0, 0], values=v, nodata_value=-1)
        self.assertEqual(grid.nodata, -1.0)
        ret = grid.set_nodata_value(np.nan)
        self.assertEqual(ret, grid)
        self.assertTrue(np.isnan(grid.nodata))
        self.assertEqual(np.sum(np.isnan(grid[:,:])), 4)
        self.assertEqual(np.sum(grid[:,:] == -1.0), 0)

#class TestInterpolation(unittest.TestCase):
#
#    def test_idw(self):
#        # Test inverse weighted distance interpolation
#        Xm, Ym = np.meshgrid(np.arange(10), np.arange(10))
#        G = np.c_[Xm.flat, Ym.flat]
#        obs = np.array(((2,4,5),(6,3,4),(3,9,9)))
#        D = raster.interpolation.interp_idw(G, obs, wfunc=lambda d: 1.0/(d+1e-16)**0.5)
#        self.assertEqual(D, np.array([ 5.77099225,  5.73080888,  5.68934588,
#                            5.64656361,  5.60263926, 5.56391335,  5.54543646,
#                            5.55810434,  5.59468599,  5.64001991, 5.7666018 ,
#                            5.71712039,  5.66733266,  5.61528705,  5.55254697,
#                            5.48200769,  5.44194711,  5.4718629 ,  5.54176929,
#                            5.61255561, 5.76627041,  5.70140932,  5.64212525,
#                            5.59011147,  5.50963694, 5.36912109,  5.24490032,
#                            5.34965772,  5.49388958,  5.59812945, 5.7765376 ,
#                            5.67823283,  5.58875282,  5.57796102,  5.51352082,
#                            5.30130025,  4.00000002,  5.2696964 ,  5.49251053,
#                            5.61373077, 5.81944097,  5.68356449,  5.00000001,
#                            5.61034065,  5.60629104, 5.47740263,  5.34297441,
#                            5.43952794,  5.57499849,  5.6693856 , 5.91915228,
#                            5.83714501,  5.76171956,  5.78893204,  5.78328667,
#                            5.71759216,  5.65687878,  5.66124481,  5.70678887,
#                            5.75517987, 6.06116315,  6.0515271 ,  6.04980015,
#                            6.05331942,  6.02324779, 5.95580534,  5.88858853,
#                            5.8514621 ,  5.84418368,  5.85242989, 6.21638838,
#                            6.27774217,  6.35281161,  6.38711869,  6.31999906,
#                            6.19926611,  6.09004119,  6.01415787,  5.96971255,
#                            5.94708184, 6.35785785,  6.4954344 ,  6.70982294,
#                            6.88076712,  6.68090612, 6.43193841,  6.26024319,
#                            6.14772366,  6.07568133,  6.03068086, 6.45468497,
#                            6.63857194,  6.9936249 ,  8.99999996,  6.97005635,
#                            6.58706897,  6.37686191,  6.24425398,  6.15677403,
#                            6.09829941]))

def peaks(n=49):
    """ 2d peaks function of MATLAB logo fame. """
    X, Y = np.meshgrid(np.linspace(-3, 3, n), np.linspace(-3, 3, n))
    return 3.0 * (1-X)**2 * np.exp(-X**2 - (Y+1)**2) \
            - 10.0 * (X/5.0 - X**3 - Y**5) * np.exp(-X**2 - Y**2) \
            - 1.0/3.0 * np.exp(-(X+1)**2 - Y**2)

if __name__ == "__main__":
    unittest.main()

