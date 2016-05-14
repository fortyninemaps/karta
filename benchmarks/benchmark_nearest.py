import numpy as np
import karta

def setup():
    t = np.linspace(0, 2*np.pi, 1000)
    coast = karta.Line(zip(3*np.sin(4*t)*np.cos(5*t), 2*np.cos(3*t)), crs=karta.crs.LonLatWGS84)

    x = np.linspace(-2, 2, 10)
    trench = karta.Line(zip(x, 0.1*x**2+3), crs=karta.crs.LonLatWGS84)
    return coast, trench

def benchmark(coast, trench):
    mindist = 1e9
    minpt = None

    for pt in coast:
        nearest_trench_pt = trench.nearest_on_boundary(pt)
        d = pt.distance(nearest_trench_pt)
        if d < mindist:
            mindist = d
            minpt = (nearest_trench_pt, pt)
    return minpt, mindist

args = setup()
res = benchmark(*args)
print(res)

