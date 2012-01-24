
import sys
sys.path.append("../")

from numpy import *
import matplotlib.pyplot as plt
import raster.flow as flow
import raster.raster as raster
import pdb, time

D = raster.witch_of_agnesi(nx=75, ny=75, a=10)

R, S = flow.dem_flow(D)

plt.imshow(S)









#upstream_area = flow.dinfinity(D)

#plt.imshow(upstream_area, cmap="gray")
#plt.ylim([0,75])
#plt.xlim([0,75])

#print upstream_area/pi

plt.show()

