
import sys
sys.path.append("../")

from numpy import *
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
import raster.flow as flow
import raster.raster as raster
import pdb, time

t0 = time.time()

D = raster.witch_of_agnesi(nx=75, ny=75, a=10)

R, S = flow.dem_flow(D)

#P = flow.build_contribution_matrix(R)

#UA = spsolve(P, ones(R.size))

#UA = flow.get_upslope_area(R)

print time.time()-t0

#pdb.set_trace()
#plt.imshow(UA, cmap="gray")

plt.show()

