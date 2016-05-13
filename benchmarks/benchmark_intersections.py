import numpy as np
from karta import Point, Multipoint, Line, Polygon
import cProfile, pstats, io

# Make a large polygon, then iterate through segments and check for intersections with a line

def benchmark():
    θ = np.linspace(0, 2*np.pi, 361)[:-1]
    r = np.sin(θ*20) + 1.5
    x = np.cos(θ) * r
    y = np.sin(θ) * r

    polygon = Polygon(zip(x, y))
    line = Line([(-2,-3), (0, 3)])

    for seg in polygon.segments:
        x = seg.intersections(line)

    bbox = Polygon([(-1, -1), (-1, 1), (1, 1), (1, -1)])

    points = []
    for pt in polygon:
        if bbox.contains(pt):
            points.append(pt)

    multipoint = Multipoint(points)
    hull = multipoint.convex_hull()

    print("bbox contains", len(multipoint), "out of", len(polygon), "vertices")
    return polygon, line, bbox, multipoint, hull

prof = cProfile.Profile()
prof.enable()

for i in range(35):
    polygon, line, bbox, multipoint, hull = benchmark()

prof.disable()
string = io.StringIO()
ps = pstats.Stats(prof, stream=string).sort_stats("tottime")
ps.print_stats(15)
print(string.getvalue())

PLOT = False
if PLOT:
    import matplotlib.pyplot as plt
    plt.plot(*polygon.coordinates)
    plt.plot(*line.coordinates)
    plt.plot(*bbox.coordinates)
    plt.plot(*multipoint.coordinates, marker=".", linewidth=0)
    plt.plot(*hull.coordinates)

    plt.show()

