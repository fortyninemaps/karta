""" Defines basic geodetic operations on a planes and spheres. """

#
# Functions are defined with an underscored version that works on scalars and a
# public version that dispatches between scalars and vectors. This approach is
# used rather than vectorization to make the C translation more straighforward.
#

from __future__ import division
import numpy as np
from math import sqrt, sin, cos, tan, asin, acos, atan, atan2, pi
import warnings

# ---------------------------------
# Functions for planar geometry
# ---------------------------------

def _plane_distance(x1, y1, x2, y2):
    return sqrt((x2 - x1)**2 + (y2 - y1)**2)

def plane_distance(xs1, ys1, xs2, ys2):
    if hasattr(xs1, "__len__"):
        return np.array([_plane_distance(x1, y1, x2, y2)
                        for x1, y1, x2, y2 in zip(xs1, ys1, xs2, ys2)])
    else:
        return _plane_distance(xs1, ys1, xs2, ys2)

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


# ---------------------------------
# Functions for spherical geodesy
# ---------------------------------

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

# ---------------------------------
# Functions for ellipsoidal geodesy
# ---------------------------------

def solve_astroid(a, f, lambda12, phi1, phi2):
    """ Used to provide an initial guess to the inverse problem in the case of
    nearly antipodal points.

    Parameters
    ----------
    a:          equatorial radius
    f:          flattening
    lambda12:   difference in longitudes (radians)
    phi1:       first latitude (radians)
    phi2:       second latitude (radians)

    Returns
    -------
    alpha1:     estimated forward azimuth at first point

    see Karney (2013) J. Geod. for details
    """
    beta1 = atan((1-f) * tan(phi1))
    beta2 = atan((1-f) * tan(phi2))
    delta = f*a*pi*cos(beta1)**2
    x = (lambda12-pi) * (a*cos(beta1)) / delta
    y = (beta2 + beta1) * a / delta
    mu = fzero_brent(1e-3, pi*a, lambda mu: mu**4 + 2*mu**3 + (1-x**2-y**2)*mu**2 - 2*y**2*mu - y**2, 1e-12)
    alpha1 = atan2(-x / (1+mu), y/mu)
    return alpha1

def solve_vicenty(a, f, lambda12, phi1, phi2):
    """ Used to provide an initial guess to the inverse problem by solving the
    corresponding problem on a sphere.

    Parameters
    ----------
    a:          equatorial radius
    f:          flattening
    lambda12:   difference in longitudes (radians)
    phi1:       first latitude (radians)
    phi2:       second latitude (radians)

    Returns
    -------
    alpha1:     forward azimuth at first point
    alpha2:     forward azimuth at second point
    s12:        distance between points

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
        alpha2 = acos(sqrt(cos(alpha1)**2*cos(beta1)**2 + (cos(beta2)**2 - cos(beta1)**2)) / cos(beta2))
    except ValueError:
        alpha2 = asin(sin(alpha0) / cos(beta2))     # Less accurate?
    sigma2 = atan2(sin(beta2), cos(alpha2)*cos(beta2))
    omega2 = atan2(sin(alpha0)*sin(sigma2), cos(sigma2))
    return alpha2, sigma2, omega2

def _normalize_longitude(x):
    """ Return longitude in the range [-180, 180). """
    return (x+180) % 360 - 180

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

    x2 = _normalize_longitude(x2-x1)
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
    a:          equatorial radius
    b:          polar radius
    x:          longitude at start
    y:          latitude at start
    azimuth:    direction travelled from point
    distnce:    distance travelled

    Returns
    -------
    x2:         longitude at destination
    y2:         latitude at destination
    back_az:    back azimuth from destination

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
    a:          equatorial radius
    b:          polar radius
    x1:         first longitude
    y1:         first latitude
    x2:         second longitude
    y2:         second latitude

    Returns
    -------
    az:         forward azimuth from first point
    back_az:    back azimuth from second point
    distance:   distance between points

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
            alpha1, _, _ = solve_vicenty(a, f, lambda12, phi1, phi2)
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


###### Utility functions ######
def unroll_angle(alpha):
    if hasattr(alpha, "__len__"):
        alpha = np.array([unroll_angle(a) for a in alpha])
    else:
        while alpha < 0:
            alpha += 2*pi
        while alpha >= 2*pi:
            alpha -= 2*pi
    return alpha

def _degrees(r):
    return r*180/pi

def fzero_brent(a, b, f, tol):
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

    raise ValueError("maximum iterations exceeded (%f, %f)" % (b, s))
