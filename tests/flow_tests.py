
import sys
sys.path.append("../..")

from numpy import *
import matplotlib.pyplot as plt
import geo_tools.raster.flow as flow

X, Y = meshgrid(range(100), range(100))

D = 1.1 * X + 0.8 * Y

upstream_area = flow.dinfinity(D)

plt.imshow(upstream_area)

#plt.imshow(D)

plt.show()

