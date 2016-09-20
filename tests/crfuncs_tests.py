import unittest
import numpy as np
from karta.raster import crfuncs

def witch_of_agnesi(nx=100, ny=100, a=4.0):
    """ Return a raster field defined by the equation Z = 8a^3 / (d^2 + 2a^2)
    where d is the distance from the center.

    Parameters
    ----------
    nx, ny : int
        raster size
    a : float
        magnitude

    Returns
    -------
    ndarray
    """
    xc = int(np.floor(nx / 2.0))
    yc = int(np.floor(ny / 2.0))
    X, Y = np.meshgrid(range(nx), range(ny))
    D = np.sqrt( (X-xc)**2 + (Y-yc)**2 )

    return (8.0 * a**3) / (D**2 + 4 * a**2)

class TestCrfuncs(unittest.TestCase):

    def test_fillarray_double(self):

        n = 50
        arr = np.zeros([n, 2*n], dtype=np.float64)
        J, I = np.meshgrid(np.arange(0, 2*n), np.arange(0, n))
        I = I.ravel().astype(np.int32)
        J = J.ravel().astype(np.int32)
        Zorig = witch_of_agnesi(2*n, n, a=7).astype(np.float64)
        Z = Zorig.ravel().copy()

        I = np.r_[I[:2232], I[2234:]]
        J = np.r_[J[:2232], J[2234:]]
        Z = np.r_[Z[:2232], Z[2234:]]

        err = crfuncs.fillarray_double(arr, I, J, Z, -999.0)
        if err != 0:
            self.fail("fillarray_double returned nonzero value")

        self.assertEqual(arr[22, 32], -999.0)
        self.assertEqual(np.sum(np.abs(Zorig[arr!=-999] - arr[arr!=-999])), 0.0)

if __name__ == "__main__":
    unittest.main()
