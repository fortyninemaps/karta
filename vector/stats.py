"""
Variogram estimation and modelling functions
"""

import numpy as np
from scipy.spatial.distance import pdist
import random

def pvariance(A):
    """ For data `A`, return the pairwise variance (squared differences). This
    will need to be Cythonized. """
    n = len(A)
    V = np.empty(sum(range(n)))
    for i in range(len(A)):
        for j in range(i):
            V[i+j] = (A[i] - A[j])**2
    return V

def estimate_vario(points, npoints=150, max_dist=None, interval=None):
    """ Given a Multipoint input, estimate a spatial variogram. By default, a
    variogram is quickly estimated. If npoints is None, then the full variogram
    is calculated.

    Parameters:
    -----------
    points : a Multipoint instance
    npoints : number of values use for variogram appromixation. If npoints is
              None, the full variogram is calculated.
    max_dist : the maximum lag
    interval : the lag interval

    Returns:
    --------
    lags : ndarray
    variogram : ndarray
    """
    if npoints is None:
        npoints = len(points.vertices)

    verts = random.sample(points.vertices, npoints)
    dist = pdist(verts)
    vari = pvariance(verts)     # Is this in the right order?

    if max_dist is None:
        max_dist = dist.max()

    if interval is None:
        interval = max_dist / 10.0

    lags = np.arange(0, max_dist+interval, interval)
    sigma_variance = np.empty_like(lags)

    for lag in lags:
        sigma_variance = np.mean(vari[dist<=lag])
    return lags, sigma_variance

def sph_func(a, h):
    """ Spherical model function for range *a* and lag *h*. """
    if h < a:
        val = 1.5*h/a - 0.5 * (h/a)**3
    else:
        val = 1
    return val

def gau_func(a, h):
    """ Gaussian model function for range *a* and lag *h*. """
    val = 1.0 - np.exp(-(h**2 / a**2))
    return val
