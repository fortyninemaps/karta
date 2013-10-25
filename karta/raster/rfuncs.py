import numpy as np

def diffrad(a, b):
    """ Return the difference in radians between two angles """
    pi = np.pi
    return (a + pi - b) % (2*pi) - pi

