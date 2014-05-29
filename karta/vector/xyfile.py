""" Convenience functions for XY(Z) tables. """

from functools import reduce
import numpy as np

def distance_xy(A):
    """ Calculate the cumulative distance at each coordinate pair in *A*. """
    d_ = np.sqrt(np.sum((A[1:,:] - A[:-1,:])**2, axis=1))
    d = [0]
    for a in d_:
        d.append(d[-1]+a)
    return np.array(d)

contains = lambda s, S: s in S

def read_xy(fnm, delimiter='', header_rows=0):
    """ Load a flowline file and return a size-2 array of coordinates. """
    notblank = lambda s: len(s.strip()) > 0
    with open(fnm) as f:
        lines = filter(notblank, f.readlines()[header_rows:])

    if delimiter == '':
        for delimiter in (',', '\t', ' '):
            if contains(delimiter, lines[-1]):
                break
            delimiter = None

    data = [[float(num) for num in line.strip().split(delimiter)]
                for line in lines]
    return np.array(data)

loadxy = read_xy
load_xy = read_xy

def write_xy(dat, fnm, delimiter=' ', header=None):
    """ Write table data in array *dat* as a delimited ASCII text file. """
    cat = lambda a,b: str(a) + delimiter + str(b)
    with open(fnm, 'w') as f:
        if header is not None:
            f.write(header+'\n')
        for line in dat:
            s = reduce(cat, line) + '\n'
            f.write(s)
    return

def xyz2array_reg(X, Y, Z):
    """ Return an array from X,Y,Z vectors, assuming gridding is
    regular. """
    xmin = min(X)
    ymin = min(Y)

    nx = sum([y==ymin for y in Y])
    ny = sum([x==xmin for x in X])

    XYZ = [(x,y,z) for x,y,z in zip(X,Y,Z)]
    XYZ.sort(key=lambda a: a[1])
    XYZ.sort(key=lambda a: a[0])
    Zs = [a[2] for a in XYZ]

    A = np.zeros([ny, nx])
    for i in range(ny):
        A[i,:] = Zs[i*nx:(i+1)*nx]

    return A


def array2xyz(A, X, Y):
    """ There are a few occasions when an XYZ list is the proper data
    format. This function converts from an array *A* with coordinates
    defined by *X* and *Y* to a list of (x,y,z) tuples.
    """
    xyz = []
    m, n = A.shape
    for j in range(n):
        for i in range (m):
            xyz.append( (X[i], Y[j], A[i,j]) )
    return xyz

