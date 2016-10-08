""" Test constructing geometry instances """

import unittest
import numpy as np
from karta import Point, Line, Polygon, Multipoint, Multiline, Multipolygon
from karta.vector.geometry import multipart_from_singleparts, merge_multiparts
from karta.crs import LonLatWGS84, WebMercator

class TestSinglepartGeometry(unittest.TestCase):

    def test_line_from_points(self):
        x = range(-5, 5)
        y = [x_**2 for x_ in x]
        vertices = list(zip(x, y))
        data = [x_*y_ for (x_, y_) in vertices]

        pts = [Point((x, y), properties={"p":i}, crs=LonLatWGS84)
                for i,((x,y),d) in enumerate(zip(vertices, data))]

        line = Line([pt.vertex for pt in pts], crs=LonLatWGS84)
        ans = Line(vertices, crs=LonLatWGS84)
        self.assertEqual(line, ans)
        return

    def test_polygon_from_points(self):
        x = range(-5, 5)
        y = [x_**2 for x_ in x]
        vertices = list(zip(x, y))
        data = [x_*y_ for (x_, y_) in vertices]

        pts = [Point((x, y), properties={"p":i}, crs=LonLatWGS84)
                for i,((x,y),d) in enumerate(zip(vertices, data))]
        poly = Polygon([pt.vertex for pt in pts], crs=LonLatWGS84)
        ans = Polygon(vertices, crs=LonLatWGS84)
        self.assertEqual(poly, ans)
        return

class TestMultipartGeometry(unittest.TestCase):

    def test_multipoint_datadict(self):
        # create a line
        vertices = [(2.0, 9.0, 9.0), (4.0, 1.0, 9.0), (4.0, 1.0, 5.0),
                    (2.0, 8.0, 0.0), (9.0, 8.0, 4.0), (1.0, 4.0, 6.0),
                    (7.0, 3.0, 4.0), (2.0, 5.0, 3.0), (1.0, 6.0, 6.0),
                    (8.0, 1.0, 0.0), (5.0, 5.0, 1.0), (4.0, 5.0, 7.0),
                    (3.0, 3.0, 5.0), (9.0, 0.0, 9.0), (6.0, 3.0, 8.0),
                    (4.0, 5.0, 7.0), (9.0, 9.0, 4.0), (1.0, 4.0, 7.0),
                    (1.0, 7.0, 8.0), (9.0, 1.0, 6.0)]

        data0 = [99.0, 2.0, 60.0, 75.0, 71.0, 34.0, 1.0, 49.0, 4.0, 36.0, 47.0,
                 58.0, 65.0, 72.0, 4.0, 27.0, 52.0, 37.0, 95.0, 17.0]

        data1 = [54.0, 40.0, 77.0, 18.0, 84.0, 91.0, 61.0, 92.0, 19.0, 42.0,
                 50.0, 25.0, 11.0, 80.0, 59.0, 56.0, 32.0, 8.0, 88.0, 76.0]

        Multipoint(vertices, data={'d0':data0, 'd1':data1})
        return

    def test_multipoint_from_points(self):
        x = range(-5, 5)
        y = [x_**2 for x_ in x]
        vertices = list(zip(x, y))
        data = [x_*y_ for (x_, y_) in vertices]

        pts = [Point((x, y), properties={"p":i}, crs=LonLatWGS84)
                for i,((x,y),d) in enumerate(zip(vertices, data))]

        mp = multipart_from_singleparts(pts)
        ans = Multipoint(vertices, data={"p":range(len(pts))}, crs=LonLatWGS84)
        self.assertEqual(mp, ans)
        return

    def test_multiline(self):
        vertices = []
        data = []
        for i in range(5):
            sub = []
            for j in range(5):
                sub.append((2*j+i, -1.5*j+2*i))
            vertices.append(sub)
            data.append(i*j)

        Multiline(vertices, data=data)
        return

    def test_multiline_from_lines_constructor(self):
        lines = []
        for i in range(5):
            sub = []
            for j in range(5):
                sub.append((2*j+i, -1.5*j+2*i))
            lines.append(Line(sub, crs=LonLatWGS84))

        g = Multiline(lines)
        for l, l_ in zip(g, lines):
            self.assertTrue(np.all(l.vertices == l_.vertices))
        self.assertEqual(g.crs, LonLatWGS84)
        return

    def test_multiline_from_lines(self):
        lines = []
        for i in range(5):
            sub = []
            for j in range(5):
                sub.append((2*j+i, -1.5*j+2*i))
            lines.append(Line(sub, properties={"d": i*j}, crs=LonLatWGS84))

        g = multipart_from_singleparts(lines)
        self.assertEqual(g.d["d"], [0, 4, 8, 12, 16])
        self.assertEqual(g.crs, LonLatWGS84)
        return

    def test_multipolygon(self):
        vertices = []
        data = []
        for i in range(5):
            sub = []
            for j in range(5):
                sub.append((2*j+i, -1.5*j+2*i))
            vertices.append([sub])
            data.append(i*j)

        Multipolygon(vertices, data=data)
        return

    def test_multipolygon_from_polygons_constructor(self):
        polys = []
        for i in range(5):
            sub = []
            for j in range(5):
                sub.append((2*j+i, -1.5*j+2*i))
            polys.append(Polygon(sub, crs=LonLatWGS84))

        g = Multipolygon(polys)
        for p, p_ in zip(g, polys):
            self.assertEqual(p.vertices, p_.vertices)
        self.assertEqual(g.crs, LonLatWGS84)
        return

    def test_multipolygon_from_polygons(self):
        polys = []
        for i in range(5):
            sub = []
            for j in range(5):
                sub.append((2*j+i, -1.5*j+2*i))
            polys.append(Polygon(sub, properties={"d": i*j}, crs=LonLatWGS84))

        g = multipart_from_singleparts(polys)
        self.assertEqual(g.d["d"], [0, 4, 8, 12, 16])
        self.assertEqual(g.crs, LonLatWGS84)
        return

    def test_empty_multpoint(self):
        mp = Multipoint([])
        self.assertEqual(len(mp), 0)
        return

    def test_empty_multiline(self):
        ml = Multiline([])
        self.assertEqual(len(ml), 0)
        return

    def test_empty_multipolygon(self):
        mp = Multipolygon([])
        self.assertEqual(len(mp), 0)
        return

    def test_multipoint_zip_init(self):
        x = range(-10, 10)
        y = [_x**2 for _x in x]
        Line(zip(x, y))
        return

    def test_merge_multipoints(self):
        mp1 = Multipoint(zip(np.arange(0, 5), np.arange(3, 6)),
                         properties={"a": 1, "b": 2},
                         data={"A": np.arange(5), "B": np.arange(10, 15)})
        mp2 = Multipoint(zip(np.arange(0, 5)+1, np.arange(3, 6)-2),
                         properties={"a": 1, "c": 3},
                         data={"A": np.arange(5), "C": np.arange(15, 20)})
        mp = merge_multiparts(mp1, mp2)
        self.assertTrue(len(mp), 10)
        self.assertEqual(set(mp.data.fields), set(["A"]))
        self.assertEqual(set(mp.properties.keys()), set(["a", "b", "c"]))
        self.assertEqual(mp.properties["a"], 1)

    def test_merge_multipoints_different_crs(self):
        mp1 = Multipoint(zip(np.arange(0, 5), np.arange(3, 6)),
                         properties={"a": 1, "b": 2},
                         data={"A": np.arange(5), "B": np.arange(10, 15)},
                         crs=LonLatWGS84)
        mp2 = Multipoint(zip(np.arange(0, 5)+1, np.arange(3, 6)-2),
                         properties={"a": 1, "c": 3},
                         data={"A": np.arange(5), "C": np.arange(15, 20)},
                         crs=WebMercator)
        mp = merge_multiparts(mp1, mp2, crs=WebMercator)
        self.assertTrue(len(mp), 10)
        self.assertEqual(set(mp.data.fields), set(["A"]))
        self.assertEqual(set(mp.properties.keys()), set(["a", "b", "c"]))
        self.assertEqual(mp.properties["a"], 1)
        self.assertEqual(mp.crs, WebMercator)
        self.assertTrue(np.allclose(mp.coordinates,
                np.array([[0.0,  1.11319491e+05, 2.22638982e+05, 1.0, 2.0, 3.0],
                          [3.34111171e+05, 4.45640110e+05, 5.57305257e+05, 1.0, 2.0, 3.0]])))

    def test_merge_multilines(self):
        x = np.arange(5)
        y = np.arange(5, 10)
        mp1 = Multiline([list(zip(x, y)), list(zip(-x, -y)), list(zip(x, -y))],
                         properties={"a": 1, "b": 2},
                         data={"A": np.arange(5), "B": np.arange(10, 15)})
        mp2 = Multiline([list(zip(x, 2*y)), list(zip(-2*x, y)), list(zip(1.5*x, y))],
                         properties={"a": 1, "c": 3},
                         data={"A": np.arange(5), "C": np.arange(15, 20)})
        mp = merge_multiparts(mp1, mp2)
        self.assertTrue(len(mp), 10)
        self.assertEqual(set(mp.data.fields), set(["A"]))
        self.assertEqual(set(mp.properties.keys()), set(["a", "b", "c"]))
        self.assertEqual(mp.properties["a"], 1)


if __name__ == "__main__":
    unittest.main()
