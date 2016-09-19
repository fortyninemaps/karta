""" Defines basic geodetic operations on a spheres and ellipsoids. """

#
# Functions are defined with an underscored version that works on scalars and a
# public version that dispatches between scalars and vectors. This approach is
# used rather than vectorization to make the C translation more straighforward.
#

from __future__ import division
import numpy as np
from math import sqrt, sin, cos, tan, asin, acos, atan, atan2, atanh, pi
from .errors import NoIntersection
import warnings

# ---------------------------------
# Miscellaneous functions
# ---------------------------------

def recurse_iterables(f, *args, **kwargs):
    def func(*args, **kwargs):
        if hasattr(args[0], "__iter__"):
            return np.array([f(*argsset, **kwargs) for argset in zip(*args)])
        else:
            return f(*args, **kwargs)
    return func

def cross(u, v):
    w = (u[1]*v[2] - u[2]*v[1],
         u[2]*v[0] - u[0]*v[2],
         u[0]*v[1] - u[1]*v[0])
    return w

def sph2cart(lon, lat):
    theta = 90-lat
    x = sin(pi*theta/180.0)*cos(pi*lon/180.0)
    y = sin(pi*theta/180.0)*sin(pi*lon/180.0)
    z = cos(pi*theta/180.0)
    return x, y, z

def cart2sph(x, y, z):
    if abs(z) > 1e-4:
        theta = atan(sqrt(x**2+y**2)/z)
    else:
        theta = acos(z / sqrt(x**2+y**2+z**2))
    if abs(x) > 1e-4:
        lon = atan2(y, x)
    else:
        lon = asin(y/sqrt(x**2+y**2))
    lat = 0.5*pi - theta
    return lon*180/pi, lat*180/pi

@recurse_iterables
def unroll_rad(rad):
    """ Return *rad* in the range [0, 2pi) """
    return rad % (2*pi)

@recurse_iterables
def reduce_rad(rad):
    """ Return *rad* in the range [-pi, pi) """
    return (rad+pi) % (2*pi) - pi

def unroll_deg(deg):
    """ Return *deg* in the range [0, 360) """
    return deg % 360

def reduce_deg(deg):
    """ Return *deg* in the range [-180, 180) """
    return (deg+180) % 360 - 180

def _degrees(r):
    return r*180.0/pi

def _radians(d):
    return d/180.0*pi

# ---------------------------------
# Functions for planar geometry
# ---------------------------------

@recurse_iterables
def plane_distance(x1, y1, x2, y2):
    return sqrt((x2 - x1)**2 + (y2 - y1)**2)

@recurse_iterables
def plane_azimuth(x1, y1, x2, y2):
    """ Cartesian azimuth between points """
    dx = x2 - x1
    dy = y2 - y1
    return atan2(dx, dy)

# ---------------------------------
# Functions for spherical geodesy
# ---------------------------------

@recurse_iterables
def sphere_distance(lons1, lats1, lons2, lats2, radius=1.0):
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

@recurse_iterables
def sphere_azimuth(lon1, lat1, lon2, lat2):
    dlon = lon2 - lon1
    return atan2(sin(_radians(dlon)), cos(_radians(lat1)) * tan(_radians(lat2)) - sin(_radians(lat1)) * cos(_radians(dlon)))

def spherical_area(r, x1, y1, x2, y2):
    """ Area between a geodesic and the equator on a sphere """
    if x2 < x1:
        reverse = -1
    else:
        reverse = 1
    _, x1, y1, x2, y2 = _canonical_configuration(x1, y1, x2, y2)
    phi1 = _radians(y1)
    phi2 = _radians(y2)
    lambda12 = _radians(x2-x1)
    alpha1, alpha2, _ = solve_vincenty(r, 0, lambda12, phi1, phi2)
    return reverse * r**2 * (alpha2-alpha1)

def isbetween_circular(x, x0, x1):
    """ Tests whether *x* is between *x0* and *x1* on a circle [-180, 180) """
    if unroll_deg(x1-x0) > 180:
        x0, x1 = x1, x0
    x = reduce_deg(x-x0)
    x1 = reduce_deg(x1-x0)
    return 0 <= x <= x1

def eulerpole(pt1, pt2):
    """ Return the Euler pole of a geodesic passing through two points. """
    v1 = sph2cart(*pt1)
    v2 = sph2cart(*pt2)
    return cross(v1, v2)

def check_in_segment_range(pt, segment):
    """ Test whether *pt* is within the lon/lat range of *segment* """
    lon_1 = unroll_deg(segment[0][0])
    lon_2 = unroll_deg(segment[1][0])

    # check longitude difference. if it's very small and latitude difference is
    # larger, go to startegy 2
    dlon = abs(reduce_deg(lon_2-lon_1))
    if (dlon < 1e-8) and (abs(segment[0][1]-segment[1][1]) > dlon):
        # strategy 2
        return (abs((pt[0]+360)%360 - lon_1) < 180) and \
               (segment[0][1] <= pt[1] <= segment[1][1])
    else:
        # strategy 1
        return isbetween_circular(pt[0], lon_1, lon_2)

def intersection_spherical(segment1, segment2):
    """ compute the intersection between two great circle segments on a sphere

    Parameters
    ----------
    segment1, segment2 : tuple of tuples
        Tuples of the form ((x1, y1), (x2, y2)) describing line segements
    """
    gc1 = eulerpole(*segment1)
    gc2 = eulerpole(*segment2)
    n = cross(gc1, gc2)
    lon, lat = cart2sph(*n)
    lon_antipodal = (lon+360)%360-180

    pt = (lon, lat)
    pt_antipodal = (lon_antipodal, -lat)
    if (check_in_segment_range(pt, segment1) and
        check_in_segment_range(pt, segment2)):
        return pt
    elif (check_in_segment_range(pt_antipodal, segment1) and
          check_in_segment_range(pt_antipodal, segment2)):
        return pt_antipodal
    else:
        raise NoIntersection()

# ---------------------------------
# Functions for ellipsoidal geodesy
# ---------------------------------

def solve_astroid(a, f, lambda12, phi1, phi2):
    """ Used to provide an initial guess to the inverse problem in the case of
    nearly antipodal points.

    Parameters
    ----------
    a: float
        equatorial radius
    f: float
        flattening
    lambda12: float (radians)
        difference in longitudes
    phi1: float (radians)
        first latitude
    phi2: float (radians)
        second latitude

    Returns
    -------
    alpha1: float (radians)
        estimated forward azimuth at first point

    see Karney (2013) J. Geod. for details
    """
    beta1 = atan((1-f) * tan(phi1))
    beta2 = atan((1-f) * tan(phi2))
    delta = f*a*pi*cos(beta1)**2
    x = (lambda12-pi) * (a*cos(beta1)) / delta
    y = (beta2 + beta1) * a / delta
    mu = fzero_brent(1e-3, pi*a, lambda mu: (mu**4 + 2*mu**3 +
                                             (1-x**2-y**2)*mu**2 - 2*y**2*mu -
                                             y**2), 1e-12)
    alpha1 = atan2(-x / (1+mu), y/mu)
    return alpha1

def solve_vincenty(a, f, lambda12, phi1, phi2):
    """ Used to provide an initial guess to the inverse problem by solving the
    corresponding problem on a sphere.

    Parameters
    ----------
    a: float
        equatorial radius
    f: float
        flattening
    lambda12: flat(radians)
        difference in longitudes
    phi1: float (radians)
        first latitude
    phi2: float (radians)
        second latitude

    Returns
    -------
    alpha1: float (radians)
        forward azimuth at first point
    alpha2: float (radians)
        forward azimuth at second point
    s12: float
        distance between points

    see Karney (2013) J. Geod. for details
    """
    eccn2 = f*(2-f)
    beta1 = atan((1-f) * tan(phi1))
    beta2 = atan((1-f) * tan(phi2))
    w = sqrt(1 - eccn2 * (0.5 * (cos(beta1) + cos(beta2)))**2)
    omega12 = lambda12 / w

    z1_r = cos(beta1)*sin(beta2) - sin(beta1)*cos(beta2)*cos(omega12)
    z1_i = cos(beta2)*sin(omega12)
    z1 = sqrt(z1_r**2 + z1_i**2)
    sigma12 = atan2(z1, sin(beta1)*sin(beta2) + cos(beta1)*cos(beta2)*cos(omega12))
    z2_r = -sin(beta1)*cos(beta2) + cos(beta1)*sin(beta2)*cos(omega12)
    z2_i = cos(beta1)*sin(omega12)

    alpha1 = atan2(z1_i, z1_r)
    alpha2 = atan2(z2_i, z2_r)
    s12 = a*w*sigma12
    return alpha1, alpha2, s12

def _solve_NEA(alpha0, alpha1, beta1):
    """ See Karney (2013) """
    sigma1 = atan2(sin(beta1), cos(alpha1)*cos(beta1))
    omega1 = atan2(sin(alpha0)*sin(sigma1), cos(sigma1))
    return sigma1, omega1

def _solve_NEB(alpha0, alpha1, beta1, beta2):
    """ See Karney (2013) """
    try:
        alpha2 = acos(sqrt(cos(alpha1)**2*cos(beta1)**2 + (cos(beta2)**2 - \
                    cos(beta1)**2)) / cos(beta2))
    except ValueError:
        alpha2 = asin(sin(alpha0) / cos(beta2))     # Less accurate?
    sigma2 = atan2(sin(beta2), cos(alpha2)*cos(beta2))
    omega2 = atan2(sin(alpha0)*sin(sigma2), cos(sigma2))
    return alpha2, sigma2, omega2

def _canonical_configuration(x1, y1, x2, y2):
    """ Put coordinates into a configuration where (Karney, eqn 44)
        y1 <= 0
        y1 <= y2 <= -y1
        0 <= x2-x1 <= 180
    """
    transformation = dict(yflip=False, xflip=False, ysignswap=False)

    if abs(y1) < abs(y2):
        y1, y2 = y2, y1
        transformation["yflip"] = True

    if y1 > 0:
        y1, y2 = -y1, -y2
        transformation["ysignswap"] = True

    x2 = reduce_deg(x2-x1)
    x1 = 0.0

    if (x2 < 0) or (x2 > 180):
        x2 = -x2
        transformation["xflip"] = True

    return transformation, x1, y1, x2, y2

def ellipsoidal_forward(a, b, x, y, azimuth, distance):
    """ Compute the destination reached starting from a point and travelling
    in a specified direction.

    Parameters
    ----------
    a: float
        equatorial radius
    b: float
        polar radius
    x: float (degrees)
        longitude at start
    y: float (degrees)
        latitude at start
    azimuth: float (radians)
        direction travelled from point
    distnce: float
        distance travelled

    Returns
    -------
    x2: float (degrees)
        longitude at destination
    y2: float (degrees)
        latitude at destination
    back_az: float (degrees)
        back azimuth from destination

    Algorithm due to Karney, C.F.F. "Algorithms for geodesics", J. Geod (2013)
    """
    f = (a-b) / a

    phi1 = pi*y/180.0
    alpha1 = pi*azimuth/180.0

    # Solve triangle NEA from Karney Fig. 1
    beta1 = atan((1-f)*tan(phi1))
    _i = sqrt(cos(alpha1)**2 + (sin(alpha1)*sin(beta1))**2)
    alpha0 = atan2(sin(alpha1)*cos(beta1), _i)
    # sigma1 = atan2(sin(beta1), cos(alpha1)*cos(beta1))
    # omega1 = atan2(sin(alpha0)*sin(sigma1), cos(sigma1))
    sigma1, omega1 = _solve_NEA(alpha0, alpha1, beta1)

    # Determine sigma2
    eccn2 = (f*(2-f))
    second_eccn2 = eccn2 / (1-eccn2)
    k2 = second_eccn2*cos(alpha0)**2

    _rad = sqrt(1+k2)
    eps = (_rad - 1) / (_rad + 1)
    A1 = 1.0/(1-eps) * (1 + eps**2/4 + eps**4/64 + eps**6/256)
    C1 = [-1.0/2*eps + 3.0/16*eps**3 - 1.0/32*eps**5,
          -1.0/16*eps**2 + 1.0/32*eps**4 - 9.0/2048*eps**6,
          -1.0/48*eps**3 + 3.0/256*eps**5,
          -5.0/512*eps**4 + 3.0/512*eps**6,
          -7.0/1280*eps**5,
          -7.0/2048*eps**6]

    I1 = A1 * (sigma1 + sum(c*sin(2*(i+1)*sigma1) for i,c in enumerate(C1)))
    s1 = I1 * b
    s2 = s1 + distance
    tau2 = s2 / (b*A1)

    C1p = [eps/2 - 9.0/32*eps**3 + 205.0/1536*eps**5,
           5.0/16*eps**2 - 37.0/96*eps**4 + 1335.0/4096*eps**6,
           29.0/96*eps**3 - 75.0/128*eps**5,
           539.0/1536*eps**4 - 2391.0/2560*eps**6,
           3467.0/7680*eps**5,
           38081.0/61440*eps**6]

    sigma2 = tau2 + sum(c*sin(2*(i+1)*tau2) for i,c in enumerate(C1p))

    # Solve triangle NEB in Karney Fig. 1
    alpha2 = atan2(sin(alpha0), cos(alpha0)*cos(sigma2))
    _j = sqrt((cos(alpha0)*cos(sigma2))**2 + sin(alpha0)**2)
    beta2 = atan2(cos(alpha0)*sin(sigma2), _j)
    omega2 = atan2(sin(alpha0)*sin(sigma2), cos(sigma2))

    # Determine lambda12
    n = f / (2.0-f)
    n2 = n*n
    A3 = 1.0 - (1.0/2-n/2)*eps - (1.0/4 + n/8 - 3.0*n2/8)*eps**2 \
        - (1.0/16 + 3.0*n/16 + n2/16)*eps**3 - (3.0/64 + n/32)*eps**4 \
        - 3.0/128*eps**5

    C3 = [(1.0/4 - n/4)*eps + (1.0/8 - n2/8)*eps**2 + (3.0/64 + 3.0*n/64 - n2/64)*eps**3 \
            + (5.0/128 + n/64)*eps**4 + 3.0/128*eps**5,
          (1.0/16 - 3.0*n/32 + n2/32)*eps**2 + (3.0/64 - n/32 - 3*n2/64)*eps**3 \
            + (3.0/128 + n/128)*eps**4 + 5.0/256*eps**5,
          (5.0/192 - 3.0*n/64 + 5.0*n2/192)*eps**3 + (3.0/128 - 5.0*n/192)*eps**4 \
            + 7.0/512*eps**5,
          (7.0/512 - 7.0*n/256)*eps**4 + 7.0*eps**5/512,
          21.0*eps**5/2560]

    I3s1 = A3 * (sigma1 + sum(c*sin(2*(i+1)*sigma1) for i,c in enumerate(C3)))
    I3s2 = A3 * (sigma2 + sum(c*sin(2*(i+1)*sigma2) for i,c in enumerate(C3)))

    lambda1 = omega1 - f*sin(alpha0)*I3s1
    lambda2 = omega2 - f*sin(alpha0)*I3s2
    lambda12 = lambda2 - lambda1

    phi2 = atan(tan(beta2) / (1-f))
    x2 = x + lambda12*180.0/pi
    if x2 >= 180.0:
        x2 -= 360.0
    y2 = phi2*180.0/pi
    backaz = (alpha2+pi)*180/pi
    x2 = (x2+180) % 360 - 180
    backaz = (backaz+180) % 360 - 180
    return x2, y2, backaz

def _ellipsoidal_inverse_equatorial(a, x1, x2):
    diff = (x2-x1 + 180) % 360 - 180
    if diff < 0:
        az = -90.0
        baz = 90.0
    else:
        az = 90.0
        baz = -90.0
    s12 = 2 * pi * a * abs(x1-x2)/360.0
    return az, baz, s12

def ellipsoidal_inverse(a, b, x1, y1, x2, y2, tol=None):
    """ Compute the shortest path (geodesic) between two points.

    Parameters
    ----------
    a: float
        equatorial radius
    b: float
        polar radius
    x1: float (degrees)
        first longitude
    y1: float (degrees)
        first latitude
    x2: float (degrees)
        second longitude
    y2: float (degrees)
        second latitude

    Returns
    -------
    az: float (degrees)
        forward azimuth from first point
    back_az: float (degrees)
        back azimuth from second point
    distance: float
        distance between points

    Algorithm due to Karney, C.F.F. "Algorithms for geodesics", J. Geod (2013)
    """
    niter = 0
    maxiter = 100
    if tol is None:
        tol = 1e-12

    if y1 == y2 == 0:
        # Equatorial case
        return _ellipsoidal_inverse_equatorial(a, x1, x2)

    # Canonical configuration
    tr, x1, y1, x2, y2 = _canonical_configuration(x1, y1, x2, y2)

    phi1 = y1*pi/180.0
    phi2 = y2*pi/180.0
    lambda12 = (x2-x1)*pi/180.0
    f = (a-b) / a

    beta1 = atan((1-f)*tan(phi1))
    beta2 = atan((1-f)*tan(phi2))

    eccn2 = f*(2-f)
    second_eccn2 = eccn2 / (1-eccn2)

    if x1 == x2:
        # Meridional case 1
        alpha0 = alpha1 = alpha2 = omega1 = omega2 = 0.0

        _i = sqrt(cos(alpha1)**2 + (sin(alpha1)*sin(beta1))**2)
        alpha0 = atan2(sin(alpha1)*cos(beta1), _i)
        sigma1, _ = _solve_NEA(alpha0, alpha1, beta1)
        _, sigma2, _ = _solve_NEB(alpha0, alpha1, beta1, beta2)

        k2 = second_eccn2
        _rad = sqrt(1+k2)
        eps = (_rad - 1) / (_rad + 1)

    elif abs(lambda12 % (2*pi) - pi) < 1e-12:
        # Meridional case 2
        if y1 + y2 > 0:
            alpha0 = alpha1 = 0.0
            alpha2 = omega1 = omega2 = pi
        else:
            alpha0 = alpha1 = omega1 = omega2 = pi
            alpha2 = 0.0

        sigma1, _ = _solve_NEA(alpha0, alpha1, beta1)
        _, sigma2, _ = _solve_NEB(alpha0, alpha1, beta1, beta2)

        k2 = second_eccn2
        _rad = sqrt(1+k2)
        eps = (_rad - 1) / (_rad + 1)

    else:
        # Newton iteration

        # Guess the azimuth
        if (abs(lambda12-pi) > 0.0087) and (abs(phi1+phi2) > 0.0087):
            # not nearly antipodal
            alpha1, _, _ = solve_vincenty(a, f, lambda12, phi1, phi2)
        else:
            alpha1 = solve_astroid(a, f, lambda12, phi1, phi2)

        dlambda12 = tol + 1

        while (abs(dlambda12) > tol) and (niter != maxiter):

            # Solve triangles
            _i = sqrt(cos(alpha1)**2 + (sin(alpha1)*sin(beta1))**2)
            alpha0 = atan2(sin(alpha1)*cos(beta1), _i)
            sigma1, omega1 = _solve_NEA(alpha0, alpha1, beta1)
            alpha2, sigma2, omega2 = _solve_NEB(alpha0, alpha1, beta1, beta2)

            # Determine lambda12
            k2 = second_eccn2 * cos(alpha0)**2
            _rad = sqrt(1+k2)
            eps = (_rad - 1) / (_rad + 1)

            n = f/(2-f)
            n2 = n*n
            A3 = 1.0 - (1.0/2 - 1.0/2*n)*eps - (1.0/4 + 1.0/8*n - 3.0/8*n2)*eps**2 \
                - (1.0/16 + 3.0/16*n + 1.0/16*n2)*eps**3 - (3.0/64 + 1.0/32*n)*eps**4 \
                - 3.0/128*eps**5

            C3 = [(1.0/4 - n/4)*eps + (1.0/8 - n2/8)*eps**2 + (3.0/64 + 3.0*n/64 - n2/64)*eps**3 \
                    + (5.0/128 + n/64)*eps**4 + 3.0/128*eps**5,
                  (1.0/16 - 3.0*n/32 + n2/32)*eps**2 + (3.0/64 - n/32 - 3*n2/64)*eps**3 \
                    + (3.0/128 + n/128)*eps**4 + 5.0/256*eps**5,
                  (5.0/192 - 3.0*n/64 + 5.0*n2/192)*eps**3 + (3.0/128 - 5.0*n/192)*eps**4 \
                    + 7.0/512*eps**5,
                  (7.0/512 - 7.0*n/256)*eps**4 + 7.0*eps**5/512,
                  21.0*eps**5/2560]

            I3s1 = A3 * (sigma1 + sum(c*sin(2*(i+1)*sigma1) for i,c in enumerate(C3)))
            I3s2 = A3 * (sigma2 + sum(c*sin(2*(i+1)*sigma2) for i,c in enumerate(C3)))

            lambda1 = omega1 - f*sin(alpha0)*I3s1
            lambda2 = omega2 - f*sin(alpha0)*I3s2
            lambda12_next = lambda2 - lambda1
            dlambda12 = lambda12_next - lambda12

            if abs(dlambda12) > tol:
                # Refine alpha1
                A1 = 1.0/(1-eps) * (1 + eps**2/4 + eps**4/64 + eps**6/256)
                C1 = [-1.0/2*eps + 3.0/16*eps**3 - 1.0/32*eps**5,
                      -1.0/16*eps**2 + 1.0/32*eps**4 - 9.0/2048*eps**6,
                      -1.0/48*eps**3 + 3.0/256*eps**5,
                      -5.0/512*eps**4 + 3.0/512*eps**6,
                      -7.0/1280*eps**5,
                      -7.0/2048*eps**6]

                I1s1 = A1 * (sigma1 + sum(c*sin(2*(i+1)*sigma1) for i,c in enumerate(C1)))
                I1s2 = A1 * (sigma2 + sum(c*sin(2*(i+1)*sigma2) for i,c in enumerate(C1)))

                A2 = (1-eps) * (1 + 1.0/4*eps**2 + 9.0/64*eps**4 + 25.0/256*eps**6)
                C2 = [1.0/2*eps + 1.0/16*eps**3 + 1.0/32*eps**5,
                      3.0/16*eps**2 + 1.0/32*eps**4 + 35.0/2048*eps**6,
                      5.0/48*eps**3 + 5.0/256*eps**5,
                      35.0/512*eps**4 + 7.0/512*eps**6,
                      63.0/1280*eps**5,
                      77.0/2048*eps**6]

                I2s1 = A2 * (sigma1 + sum(c*sin(2*(i+1)*sigma1) for i,c in enumerate(C2)))
                I2s2 = A2 * (sigma2 + sum(c*sin(2*(i+1)*sigma2) for i,c in enumerate(C2)))

                Js1 = I1s1 - I2s1
                Js2 = I1s2 - I2s2

                m12 = b * (sqrt(1 + k2*sin(sigma2)**2) * cos(sigma1)*sin(sigma2) \
                         - sqrt(1 + k2*sin(sigma1)**2) * sin(sigma1)*cos(sigma2) \
                         - cos(sigma1) * cos(sigma2) * (Js2-Js1))
                dlambda12_dalpha1 = m12/(a * cos(alpha2)*cos(beta2))
                dalpha1 = -dlambda12 / dlambda12_dalpha1
                alpha1 = (alpha1 + dalpha1) % (2*pi)

            niter += 1

    if niter == maxiter:
        warnings.warn(
            "Convergence failure (%f, %f) -> (%f, %f)" % (x1, y1, x2, y2),
            RuntimeWarning)

    k2 = second_eccn2 * cos(alpha0)**2
    _rad = sqrt(1+k2)
    eps = (_rad - 1) / (_rad + 1)

    # Determine s12
    A1 = 1.0/(1-eps) * (1 + eps**2/4 + eps**4/64 + eps**6/256)
    C1 = [-1.0/2*eps + 3.0/16*eps**3 - 1.0/32*eps**5,
          -1.0/16*eps**2 + 1.0/32*eps**4 - 9.0/2048*eps**6,
          -1.0/48*eps**3 + 3.0/256*eps**5,
          -5.0/512*eps**4 + 3.0/512*eps**6,
          -7.0/1280*eps**5,
          -7.0/2048*eps**6]

    I1s1 = A1 * (sigma1 + sum(c*sin(2*(i+1)*sigma1) for i,c in enumerate(C1)))
    I1s2 = A1 * (sigma2 + sum(c*sin(2*(i+1)*sigma2) for i,c in enumerate(C1)))

    s1 = I1s1*b
    s2 = I1s2*b
    s12 = s2-s1

    if tr["xflip"]:
        alpha1 = -alpha1
        alpha2 = -alpha2

    if tr["yflip"]:
        alpha1, alpha2 = pi-alpha2, pi-alpha1

    if tr["ysignswap"]:
        alpha1 = pi - alpha1
        alpha2 = pi - alpha2

    az = (alpha1*180/pi + 180) % 360 - 180
    backaz = ((alpha2+pi)*180/pi + 180) % 360 - 180
    return az, backaz, s12

def _ellipsoidal_area(a, b, lambda12, phi1, phi2, alpha1, alpha2):
    """ Area of a single quatrilateral defined by two meridians, the equator,
    and another geodesic.
    """
    f = (a-b)/a
    e2 = f*(2-f)
    ep2 = e2/(1-e2)
    e = sqrt(e2)

    # Authalic radius
    c = sqrt(a**2/2 + b**2/2*atanh(e)/e)

    beta1 = atan((1-f)*tan(phi1))
    beta2 = atan((1-f)*tan(phi2))

    _i = sqrt(cos(alpha1)**2 + (sin(alpha1)*sin(beta1))**2)
    alpha0 = atan2(sin(alpha1)*cos(beta1), _i)

    sigma1, omega1 = _solve_NEA(alpha0, alpha1, beta1)
    _, sigma2, omega2 = _solve_NEB(alpha0, alpha1, beta1, beta2)
    omega12 = omega2 - omega1

    # Bessel identity for alpha2 - alpha1
    alpha12 = 2*atan(sin(0.5*beta1+0.5*beta2)/cos(0.5*beta2-0.5*beta1) \
            * tan(0.5*omega12))
    sph_term = c**2 * alpha12

    # compute integrals for ellipsoidal correction
    k2 = ep2*cos(alpha0)**2

    C40 = (2.0/3 - ep2/15 + 4*ep2**2/105 - 8*ep2**3/315 + 64*ep2**4/3465 - 128*ep2**5/9009) \
        - (1.0/20 - ep2/35 + 2*ep2**2/105 - 16*ep2**3/1155 + 32*ep2**4/3003) * k2 \
        + (1.0/42 - ep2/63 + 8*ep2**2/693 - 90*ep2**3/9009) * k2**2 \
        - (1.0/72 - ep2/99 + 10*ep2**2/1287) * k2**3 \
        + (1.0/110 - ep2/143) * k2**4 - k2**5/156

    C41 = (1.0/180 - ep2/315 + 2*ep2**2/945 - 16*ep2**3/10395 + 32*ep2**4/27027) * k2 \
        - (1.0/252 - ep2/378 + 4*ep2**2/2079 - 40*ep2**3/27027) * k2**2 \
        + (1.0/360 - ep2/495 + 2*ep2**2/1287) * k2**3 \
        - (1.0/495 - 2*ep2/1287) * k2**4 + 5*k2**5/3276

    C42 = (1.0/2100 - ep2/3150 + 4*ep2**2/17325 - 8*ep2**3/45045) * k2**2 \
        - (1.0/1800 - ep2/2475 + 2*ep2**2/6435) * k2**3 \
        + (1.0/1925 - 2*ep2/5005) * k2**4 - k2**5/2184

    C43 = (1.0/17640 - ep2/24255 + 2*ep2**2/63063) * k2**3 \
        - (1.0/10780 - ep2/14014) * k2**4 + 5*k2**5/45864

    C44 = (1.0/124740 - ep2/162162) * k2**4 - 1*k2**5/58968

    C45 = k2**5/792792

    Cs = [C40, C41, C42, C43, C44, C45]
    I4s1 = sum(c*cos((2*i+1)*sigma1) for i,c in enumerate(Cs))
    I4s2 = sum(c*cos((2*i+1)*sigma2) for i,c in enumerate(Cs))

    S12 = sph_term + e**2*a**2 * cos(alpha0)*sin(alpha0) * (I4s2-I4s1)
    return S12


def ellipsoidal_area(a, b, x1, y1, x2, y2):
    """ Area of a single quadrilateral defined by two meridians, the equator,
    and another geodesic.

    Parameters
    ----------
    a : float
        equatorial radius
    b : float
        polar radius
    x1 : float (degrees)
        longitude of geodesic segment start
    y1 : float (degrees)
        latitude of geodesic segment start
    x2 : float (degrees)
        longitude of geodesic segment end
    y2 : float (degrees)
        latitude of geodesic segment end

    """
    if x2 < x1:
        reverse = -1
    else:
        reverse = 1
    _, x1, y1, x2, y2 = _canonical_configuration(x1, y1, x2, y2)
    phi1 = y1*pi/180.0
    phi2 = y2*pi/180.0
    lambda12 = (x2-x1)*pi/180.0

    az, baz, _ = ellipsoidal_inverse(a, b, x1, y1, x2, y2)
    alpha1 = az * pi/180
    alpha2 = (baz-pi) * pi/180
    return reverse * _ellipsoidal_area(a, b, lambda12, phi1, phi2, alpha1, alpha2)


###### Root-finding ######

def fzero_brent(a, b, f, tol):
    """ Evaluate a function to find a bracketed root

    Parameters
    ----------
    a : float
        left bracket
    b : float
        right bracket
    f : callable
        function to find zero of
    tol : float
        error tolerance

    Raises
    ------
    ValueError
        if root is not bracketed
    RuntimeError
        if maximum iterations exceeded
    """
    fa = f(a)
    fb = f(b)
    if fa == 0:
        return a
    elif fb == 0:
        return b
    elif fa * fb > 0:
        raise ValueError("root not bracketed")
    if abs(fa) < abs(fb):
        a,b = b,a
    c = a
    d = np.nan
    mflag = True

    delta = 1e-10
    niter = 0
    while niter != 1000:

        fc = f(c)
        if (fa != fc) and (fb != fc):
            # inverse quadratic interpolation
            s = a*fb*fc / ((fa-fb) * (fa-fc)) \
                    + b*fa*fc / ((fb-fa) * (fb-fc)) \
                    + c*fa*fb / ((fc-fa) * (fc-fb))
        else:
            # secant
            s = b - fb*(b-a)/(fb-fa)

        if ((not 0.25*(3*a+b) < s < b) and (not b < s < 0.25*(3*a+b))) or \
                (mflag and (abs(s-b) >= 0.5*abs(b-c))) or \
                (~mflag and (abs(s-b) >= 0.5*abs(c-d))) or \
                (mflag and (abs(b-c) < abs(delta))) or \
                (~mflag and (abs(c-d) < abs(delta))):
            # bisection
            s = 0.5*(a+b)
            mflag = True
        else:
            mflag = False

        fs = f(s)
        d = c
        c = b
        fa = f(a)
        fb = f(b)
        if fa*fs < 0:
            b = s
            fb = fs
        else:
            a = s
            fa = fs

        if abs(fa) < abs(fb):
            a,b = b,a
            fa,fb = fb,fa

        if fb == 0:
            return b
        elif f(s) == 0:
            return s
        elif (abs(b-a) < tol):
            return b
        niter += 1

    raise RuntimeError("maximum iterations exceeded (%f, %f)" % (b, s))
