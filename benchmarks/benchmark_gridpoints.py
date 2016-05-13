import karta
import numpy as np
import cProfile, pstats, io

n = 50000
np.random.seed(49)
x = np.random.rand(n)
y = np.random.rand(n)
z = (x-0.5)**2 + (y-0.2)**3 + 0.5*np.random.rand(n)

T = [0.0, 0.0, 0.01, 0.01, 0.0, 0.0]

prof = cProfile.Profile()
prof.enable()

for i in range(10):
    grid = karta.raster.gridpoints(x, y, z, T, karta.crs.Cartesian)

prof.disable()
string = io.StringIO()
sortby = "tottime"
ps = pstats.Stats(prof, stream=string).sort_stats(sortby)
ps.print_stats(15)
print(string.getvalue())
