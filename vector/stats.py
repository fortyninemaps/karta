"""
Geostatistics module

Contains variogram calculation and modelling tools
"""

import math
import numpy as np
import itertools

def calculate_vario(points, npoints=150, max_dist=None, interval=None):
    """ Given a list of Points, calculate a spatial variogram. By default,
    a variogram is quickly estimated. If npoints is None, then the full
    variogram is calculated.

    Parameters:
        points      : list of Point objects
        npoints     : number of values use for variogram appromixation
                      if npoints is None, the full variogram is calculated
        max_dist    : the maximum lag
        interval    : the lag interval

    Returns
        bins        : ndarray
        variogram   : ndarray
    """
    if isinstance(points, np.ndarray) == False:
        points = np.array(points)

    if npoints is None:
        npoints = len(points)

    # Get a list of the points that will be calculated
    ipoints = np.random.permutation(len(points))[:npoints]

    # Get a pool of points to use for sampling
    sample_points = np.random.permutation(points)

    # Fill in a datastructure for holding the raw data
    # This will be a list of lists, npoints long
    # Each sublist, corresponding to a point in ipoints,
    # will contain tuples, with every additional point's distance
    # and the variance
    V = []
    for i in ipoints:
        np.random.shuffle(sample_points)
        V.append([])
        pt1 = points[i]
        i2 = 0

        # Find all other points within max_dist of point
        for pt2 in sample_points:

            if i2 > npoints:
                # Have looked through npoints already; move on
                break

            if pt1 is pt2:
                # Skip self, so lag != 0
                pass

            else:
                d = math.sqrt((pt1.x-pt2.x)**2 + (pt1.y-pt2.y)**2)
                if (max_dist is None) or (d <= max_dist):
                    V[-1].append( (d, (pt2.z - pt1.z)**2) )

            i2 += 1

    # Initialize a set of bins for holding variogram values
    if max_dist is None:
        max_dist = max((i[0] for i in list(itertools.chain(V))))

    if interval is None:
        interval = max_dist / 10.0
    bins = np.arange(0, max_dist+interval, interval)

    # From the raw data, calculate a variogram
    variogram = []
    for b0, bf in zip(bins[:-1], bins[1:]):

        p = filter(lambda a: (a[0] >= b0) and (a[0] < bf),
                   itertools.chain(*V))
        variance = sum([i[1] for i in p]) / len(p)
        variogram.append(variance)

    return bins[:-1], np.array(variogram)


def sph_func(a, h):
    """ Spherical model function for range *a* and lag *h*. """
    if h<a:
        val = 1.5*h/a - 0.5 * (h/a)**3
    else:
        val = 1
    return val

def gau_func(a, h):
    """ Gaussian model function for range *a* and lag *h*. """
    val = 1.0 - np.exp(-(h**2 / a**2))
    return val
