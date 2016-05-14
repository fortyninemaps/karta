import timeit

res = timeit.timeit(stmt="karta.Point((2,3))",
                    setup="import karta",
                    number=100000)
print("Point: {0}".format(res))

res = timeit.timeit(stmt="karta.Line(verts)",
                    setup="import karta; import random; verts = [(random.random(), random.random()) for _ in range(500)]",
                    number=1000)
print("Line: {0}".format(res))

res = timeit.timeit(stmt="karta.Polygon(verts)",
                    setup="import karta; import random; verts = [(random.random(), random.random()) for _ in range(500)]",
                    number=1000)
print("Polygon: {0}".format(res))

res = timeit.timeit(stmt="karta.Multipoint(verts)",
                    setup="import karta; import random; verts = [(random.random(), random.random()) for _ in range(500)]",
                    number=10000)
print("Multipoint: {0}".format(res))

res = timeit.timeit(stmt="karta.Multiline(verts)",
                    setup="import karta; import random; verts = [[(random.random(), random.random()) for _ in range(100)] for _ in range(50)]",
                    number=1000)
print("Multiline: {0}".format(res))


res = timeit.timeit(stmt="karta.Multipolygon(verts)",
                    setup="import karta; import random; verts = [[[(random.random(), random.random())] for _ in range(100)] for _ in range(10)]",
                    number=1000)
print("Multipolygon: {0}".format(res))

