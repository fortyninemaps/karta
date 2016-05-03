""" Unit tests for vector functions """

import unittest
import numpy as np
import shapely.geometry

import karta.vector as vector
from karta.vector.geometry import Point, Multipoint, Line, Polygon

class TestGeoInterface(unittest.TestCase):

    def test_point(self):
        pt = Point((1,2))
        self.assertEqual(pt.geomdict, {"type":"Point", "coordinates":(1,2)})
        pt = pt.shift((2,2))
        self.assertEqual(pt.geomdict, {"type":"Point", "coordinates":(3,4)})

    def test_poly(self):
        x = np.arange(5)
        y = x**2
        poly = Polygon(list(zip(x, y)))
        self.assertEqual(poly.geomdict,
                         {"type":"Polygon",
                          "bbox":(0, 0, 4, 16),
                          "coordinates": [list(zip(x, y))]})

    def test_line(self):
        x = np.arange(10)
        y = x**2
        line = Line(zip(x, y))
        self.assertEqual(line.geomdict,
                         {"type":"LineString",
                          "bbox":(0, 0, 9, 81),
                          "coordinates": list(zip(x,y))})
        line = line[:5]
        self.assertEqual(line.geomdict,
                         {"type":"LineString",
                          "bbox":(0, 0, 4, 16),
                          "coordinates": list(zip(x[:5],y[:5]))})


    def test_point_output(self):
        p = Point((4, 2))
        sp = shapely.geometry.shape(p.geomdict)
        self.assertEqual(sp.x, p.x)
        self.assertEqual(sp.y, p.y)
        return

    def test_multipoint_output(self):
        p = Multipoint([(4, 2), (3, 5), (3, 2), (7, 3)])
        sp = shapely.geometry.shape(p.geomdict)
        x, y = p.coordinates
        self.assertEqual(x, tuple([el.x for el in sp]))
        self.assertEqual(y, tuple([el.y for el in sp]))
        return

    def test_line_output(self):
        p = Line([(4, 2), (3, 5), (3, 2), (7, 3)])
        sp = shapely.geometry.shape(p.geomdict)
        x, y = p.coordinates
        sx, sy = sp.xy
        self.assertEqual(x, tuple(sx))
        self.assertEqual(y, tuple(sy))
        return

    def test_poly_output(self):
        p = Polygon([(4, 2), (3, 5), (3, 2), (7, 3)])
        sp = shapely.geometry.shape(p.geomdict)
        self.assertEqual(p.bbox, sp.bounds)
        return

    def test_point_input(self):
        sp = shapely.geometry.Point((3,4))
        p = vector.read.from_shape(sp)
        self.assertEqual(p.x, sp.x)
        self.assertEqual(p.y, sp.y)
        return

    def test_line_input(self):
        sp = shapely.geometry.LineString([(3,4), (6,2), (2,5)])
        p = vector.read.from_shape(sp)
        x, y = p.coordinates
        sx, sy = sp.xy
        self.assertEqual(x, tuple(sx))
        self.assertEqual(y, tuple(sy))
        return

    def test_poly_input(self):
        sp = shapely.geometry.Polygon([(4, 2), (3, 5), (3, 2), (7, 3)])
        p = vector.read.from_shape(sp)
        self.assertEqual(p.bbox, sp.bounds)
        return

    def test_multipoly_input(self):
        sp1 = shapely.geometry.Polygon([(4, 2), (3, 5), (3, 2), (7, 3)])
        sp2 = shapely.geometry.Polygon([(7, 3), (9, 7), (2, 7), (2, 0)])
        smp = shapely.geometry.MultiPolygon([sp1, sp2])
        mpoly = vector.read.from_shape(smp)
        p1, p2 = mpoly
        self.assertEqual(p1.bbox, sp1.bounds)
        self.assertEqual(p2.bbox, sp2.bounds)
        return

    def test_multiline_input(self):
        sp1 = shapely.geometry.LineString([(4, 2), (3, 5), (3, 2), (7, 3)])
        sp2 = shapely.geometry.LineString([(7, 3), (9, 7), (2, 7), (2, 0)])
        smp = shapely.geometry.MultiLineString([sp1, sp2])
        p1, p2 = vector.read.from_shape(smp)
        x, y = p1.coordinates
        sx, sy = sp1.xy
        self.assertEqual(x, tuple(sx))
        self.assertEqual(y, tuple(sy))
        x, y = p2.coordinates
        sx, sy = sp2.xy
        self.assertEqual(x, tuple(sx))
        self.assertEqual(y, tuple(sy))
        return

    def test_feature_input(self):
        class Pointy(object):
            __geo_interface__ = {'type': 'Point', 'coordinates': (0.0, 0.0)}
        class Placemark(object):
            __geo_interface__ = {
                'type': 'Feature',
                'properties': {'name': 'Phoo'},
                'geometry': Pointy.__geo_interface__ }
        p = vector.read.from_shape(Placemark())
        self.assertEqual(p.properties["name"], "Phoo")
        self.assertEqual(p.vertex, (0.0, 0.0))
        return

if __name__ == "__main__":
    unittest.main()
