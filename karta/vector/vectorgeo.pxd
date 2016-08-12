""" Header file for vector geometry primitives """

cdef struct Vector2:
    double x
    double y

cdef struct Vector3:
    double x
    double y
    double z

cdef inline double dot2(Vector2, Vector2) nogil
cdef inline double dot3(Vector3, Vector3) nogil
cdef inline double cross2(Vector2, Vector2) nogil
cdef inline Vector3 cross3(Vector3, Vector3)

cdef Vector2 proj2(Vector2, Vector2)
cdef inline double dist2(Vector2, Vector2) nogil
# cdef double distsph(Vector2, Vector2)       # TODO

cdef double mind(double, double) nogil
cdef double maxd(double, double) nogil
cdef double absd(double) nogil

# cdef Vector3 euler_pole(Vector2, Vector2)   # TODO
# cdef Vector2 sph2cart(Vector3)              # TODO
# cdef Vector3 cart2sph(Vector2)              # TODO

