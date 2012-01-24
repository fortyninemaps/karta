
import sys
sys.path.append("../")

from numpy import *
import matplotlib.pyplot as plt
import raster.flow as flow
import raster.raster as raster
import pdb


#X, Y = meshgrid(range(250), range(250))

#D = 0.0 * X + 0.5 * Y
#D = 0.5 * X + 0.0 * Y
#D = 0.5 * X + 0.5 * Y
D = raster.witch_of_agnesi(nx=75, ny=75, a=10)

Dp = raster.pad(D)
c = [(i,j) for i in range(1, Dp.shape[0]-1) for j in range(1, Dp.shape[1]-1)]

ff = map(lambda a: flow.pixel_flow(Dp, a[0], a[1], d1=20.0, d2=20.0), c)

#pdb.set_trace()
R = array([i[0] for i in ff]).reshape(D.shape)
S = array([i[1] for i in ff]).reshape(D.shape)

#pdb.set_trace()

#plt.imshow(S)

#upstream_area = flow.dinfinity(D)

#plt.imshow(upstream_area, cmap="gray")
#plt.ylim([0,75])
#plt.xlim([0,75])

#print upstream_area/pi

#plt.show()

