#! /usr/bin/env python
#
#   Guppy tests
#


from guppy import *
from math import pi
import pdb

print 'simple class attributes and methods'
A = polygon([(0,1),(1,1),(1,0),(0,0.5)])
A = polygon([(0,1,0.0),(1,1,1.0),(1,0,2.0),(0,0.5,1)])
print 'rank', A.rank
print 'perimeter', A.perimeter()
print 'area', A.area()
As = A.to_shapely()
print type(As)
A.shift((2,1.5,.2))
A.print_vertices()

B = polyline([(0,1),(1,1),(1,0),(0,0.5)])
print 'length', B.length()
print 'displacement', B.displacement()
B.print_vertices()
Bs = B.to_shapely()
print type(Bs)

print
print 'walk test'
P = point((2,3))
print 'walk', P.walk(4, 3.1415)

print
print 'tighten test'
X = [1.0,2.,4.,6.,8.,10.,12.,14.,16.,18.]
Z = [2.5,3.1,2.5,3.2,3.7,3.4,3.6,2.9,2.4,2.7]
meas = zip(X,Z)
meas_t = tighten(X,Z)
print meas_t[-1]

print
print 'intersection point test'
sol = ray_intersection((2,2), (4,0), (4,4), direction=0.0)
print sol

print
print 'multipoint topology tests'
square = polygon([(0,0), (4,0), (4,4), (0,4)])
print square.contains(point((2,2)))
print square.contains(point((2,6)))
pt = point((3,3))
print square.nearest_to(pt)

endpt1 = (1,1)
endpt2 = (5,6)
pt = (0,1)
print pt_nearest(pt, endpt1, endpt2)
