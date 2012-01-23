
import sys
sys.path.append("../")

from numpy import *
import matplotlib.pyplot as plt
import raster.flow as flow
import raster.raster as raster

#X, Y = meshgrid(range(250), range(250))

#D = 0.0 * X + 0.5 * Y
#D = 0.5 * X + 0.0 * Y
#D = 0.5 * X + 0.5 * Y
D = raster.witch_of_agnesi(nx=75, ny=75)

upstream_area = flow.dinfinity(D)

plt.imshow(upstream_area)

#plt.imshow(D)

print upstream_area/pi

plt.show()

