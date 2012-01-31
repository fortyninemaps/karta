
import sys
sys.path.append("../")

from numpy import *
import matplotlib.pyplot as plt
import raster.flow as flow
import raster.raster as raster

A = raster.witch_of_agnesi(a=8.0)

Aslope = raster.slope(A)

Aaspect = raster.aspect(A)

Aflow = flow.dinfinity(A)

fig = plt.figure()
ax1 = fig.add_subplot(2,2,1)
ax2 = fig.add_subplot(2,2,2)
ax3 = fig.add_subplot(2,2,3)
ax4 = fig.add_subplot(2,2,4)

ax1.imshow(A)
ax2.imshow(Aslope)
ax3.imshow(Aaspect)

ax1.set_title("witch of agnesi")
ax2.set_title("slope")
ax3.set_title("aspect")

plt.show()
