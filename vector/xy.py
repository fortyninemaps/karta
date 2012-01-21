""" Convenience functions for XY(Z) tables. """

import numpy as np

def load_xy(fnm, delimiter=","):
    """ Load a flowline file and return a size-2 array of coordinates. """
    with open(fnm) as f:
        coords = [[float(j) for j in i.split(delimiter)] for i in f.readlines()]
    return np.array(coords)
