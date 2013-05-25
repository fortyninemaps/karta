"""
Variogram estimation and modelling functions
"""

import numpy as np
from scipy.spatial.distance import pdist
from scipy.optimize import minimize
import random

class VariogramFunction(object):
    """ Base class for variogram functions, which can be composed
    through summation. """

    def __call__(self, lags):
        return self.func(lags)

    def __add__(self, other):
        return ComposedFunction(lambda lags: self(lags) + other(lags))

    def __rmul__(self, other):
        return ComposedFunction(lambda lags: other * self(lags))

class ComposedFunction(VariogramFunction):
    """ Type returned from function algebra """
    def __init__(self, f):
        self.func = lambda h: f(h)
        return

class Sph(VariogramFunction):
    """ Spherical model function """
    def __init__(self, rng):
        self.func = lambda h: np.where(h < rng, 1.5*h/rng - 0.5 * (h/rng)**3, 1.0)
        return

class Gau(VariogramFunction):
    """ Gaussian model function """
    def __init__(self, rng):
        self.func = lambda h: 1.0 - np.exp(-(h**2 / rng**2))
        return

class Nug(VariogramFunction):
    """ Nugget effect function """
    def __init__(self):
        self.func = lambda h: 1.0
        return

def fit_model(model, p, lags, variance):
    """ Fit a variogram model to an experimental variogram.

    Parameters:
    -----------

    model : function for creating variogram estimators that takes *lags* as an
            argument, e.g.

        model = lambda p: p[0] * Nug() +  p[1] * Sph(p[2])

    p : array of arguments for *model*
    lags : variogram offsets
    variance : variogram variance corresponding to lags
    """
    lags = np.array(lags)
    obj = lambda p: sum((model(p)(lags) - variance)**2)
    res = minimize(obj, p)
    return res.x

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

    lags = np.arange(interval, max_dist+interval, interval)
    sigma_variance = np.empty_like(lags)
    for i, lag in enumerate(lags):
        band = dist <= lag
        sigma_variance[i] = np.mean(vari[band])
    return lags, sigma_variance

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

