""" Header file for vector geometry primitives """

cdef struct Vector2:
    double x
    double y

cdef struct Vector3:
    double x
    double y
    double z

cdef int fsign(double)

cdef double dot2(Vector2, Vector2) nogil
cdef double dot3(Vector3, Vector3) nogil
cdef double cross2(Vector2, Vector2) nogil
cdef Vector3 cross3(Vector3, Vector3)

cdef Vector2 proj2(Vector2, Vector2)
cdef double dist2(Vector2, Vector2) nogil
cdef double dist_sph(Vector2, Vector2) nogil
cdef double azimuth(Vector2, Vector2)
cdef double azimuth_sph(Vector2, Vector2)
cdef int bndlat_sph(Vector2, Vector2, double*, double*)

cdef Vector3 eulerpole(Vector2, Vector2)
cdef Vector3 eulerpole_cart(Vector3, Vector3)
cdef Vector3 sph2cart(Vector2)
cdef Vector2 cart2sph(Vector3)

