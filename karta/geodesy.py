""" Defines basic geodetic operations on a planes and spheres. """

#
# Functions are defined with an underscored version that works on scalars and a
# public version that dispatches between scalars and vectors. This approach is
# used rather than vectorization to make the C translation more straighforward.
#

import numpy as np
from math import sqrt, sin, cos, tan, asin, acos, atan, pi

def _plane_distance(x1, y1, x2, y2):
    return sqrt((x2 - x1)**2 + (y2 - y1)**2)

def plane_distance(xs1, ys1, xs2, ys2):
    if hasattr(xs1, "__len__"):
        return np.array([_plane_distance(x1, y1, x2, y2)
                        for x1, y1, x2, y2 in zip(xs1, ys1, xs2, ys2)])
    else:
        return _plane_distance(xs1, ys1, xs2, ys2)

def _sphere_distance(lon1, lat1, lon2, lat2, radius):
    dx = abs(lon1-lon2)
    dy = abs(lat1-lat2)

    if dx > 0.01 or dy > 0.01:
        # use spherical law of cosines
        d_ = acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(dx))

    else:
        # use haversine
        d_ = (2 * asin(sqrt(sin(dy / 2.)**2 +
                        cos(lat1) * cos(lat2) * sin(dx / 2.)**2)))
    return radius * d_

def sphere_distance(lons1, lats1, lons2, lats2, radius):
    if hasattr(lons1, "__len__"):
        return np.array([_sphere_distance(ln1, lt1, ln2, lt2, radius)
                        for ln1, lt1, ln2, lt2
                        in zip(lons1, lats1, lons2, lats2)])
    else:
        return _sphere_distance(lons1, lats1, lons2, lats2, radius)


def _plane_azimuth(x1, y1, x2, y2):
    """ Return cartesian azimuth between points (scalar func) """
    dx = x2 - x1
    dy = y2 - y1
    if dx > 0:
        if dy > 0:
            return np.arctan(dx/dy)
        elif dy < 0:
            return np.pi - np.arctan(-dx/dy)
        else:
            return 0.5*np.pi
    elif dx < 0:
        if dy > 0:
            return 2*np.pi - np.arctan(-dx/dy)
        elif dy < 0:
            return np.pi + np.arctan(dx/dy)
        else:
            return 1.5*np.pi
    else:
        if dy > 0:
            return 0.0
        else:
            return np.pi

def plane_azimuth(xs1, ys1, xs2, ys2):
    """ Return cartesian azimuth between points """
    if hasattr(xs1, "__len__"):
        return np.array([_plane_azimuth(x1, y1, x2, y2)
                        for x1, y1, x2, y2 in zip(xs1, ys1, xs2, ys2)])
    else:
        return _plane_azimuth(xs1, ys1, xs2, ys2)

def _sphere_azimuth(lon1, lat1, lon2, lat2):
    dlon = lon2 - lon1
    if (cos(lat1) * tan(lat2) - sin(lat1) * cos(dlon)) == 0:
        az = 0.5*pi
    else:
        az = atan(sin(dlon) / (cos(lat1) * tan(lat2) - sin(lat1) * cos(dlon)))
    
    if unroll_angle(dlon) <= pi:
        if lat2 < lat1:
            az = az + pi
    else:
        if lat1 < lat2:
            az = az + 2*pi
        else:
            az = az + pi
    return unroll_angle(az)

def sphere_azimuth(lons1, lats1, lons2, lats2):
    if hasattr(lons1, "__len__"):
        return np.array([_sphere_azimuth(lon1, lat1, lon2, lat2)
                        for lon1, lat1, lon2, lat2
                        in zip(lons1, lats1, lons2, lats2)])
    else:
        return _sphere_azimuth(lons1, lats1, lons2, lats2)


###### Utility functions ######
def unroll_angle(alpha):
    if hasattr(alpha, "__len__"):
        alpha_unrolled = np.array([unroll_angle(a) for a in alpha])
    else:
        while alpha < 0:
            alpha += 2*pi
        while alpha >= 2*pi:
            alpha -= 2*pi
    return alpha

