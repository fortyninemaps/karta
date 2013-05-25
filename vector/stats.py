"""
Variogram estimation and modelling functions
"""

import numpy as np
from scipy.spatial.distance import pdist
from scipy.optimize import fmin
import random

def pvariance(A):
    """ For data `A`, return the pairwise variance (squared differences). This
    will need to be Cythonized. """
    n = len(A)
    V = np.empty(sum(range(n)))
    ind = 0
    for i in range(len(A)):
        for j in range(i+1, n):
            V[ind] = (A[i] - A[j])**2
            ind += 1
    return V

def fit_model(model, p, lags, variance):
    """ Fit callable *model* with guess *p* to an experimental variogram
    given by *lags*, *variance*.

    *model* should return an estimator function that takes *lags* as an
    argument
    """
    obj = lambda p: model(p)(lags)
    res = fmin(obj, p)
    pbest = res[0]
    return

def estimate_vario(mp, npoints=150, max_dist=None, interval=None):
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
        npoints = len(mp)

    irand = random.sample(np.arange(len(mp)), npoints)
    verts = mp.get_vertices()[irand]
    z = mp.get_data()[irand]
    dist = pdist(verts)
    vari = pvariance(z)

    if max_dist is None:
        max_dist = dist.max()

    if interval is None:
        interval = max_dist / 10.0

    lags = np.arange(0, max_dist+interval, interval)
    sigma_variance = np.empty_like(lags)

    sigma_variance = [np.mean(vari[dist<=lag]) for lag in lags]
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
