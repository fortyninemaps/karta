"""
2D raster functions
"""

import numpy as np

def witch_of_agnesi(nx=100, ny=100, a=4.0):
    """ Return a raster field defined by the equation
                    8 a^3
            Z = --------------
                 d^2 + 4 a*2
    where d is the distance from the center.
    """
    xc = int(np.ceil(nx / 2.0))
    yc = int(np.ceil(ny / 2.0))
    X, Y = np.meshgrid(range(nx), range(ny))
    D = np.sqrt( (X-xc)**2 + (Y-yc)**2 )

    return (8.0 * a**3) / (D**2 + 4 * a**2)


def peaks(n=49):
    """ 2d peaks function of MATLAB logo fame. """
    X, Y = np.meshgrid(np.linspace(-3, 3, n), np.linspace(-3, 3, n))
    return 3.0 * (1-X)**2 * np.exp(-X**2 - (Y+1)**2) \
            - 10.0 * (X/5.0 - X**3 - Y**5) * np.exp(-X**2 - Y**2) \
            - 1.0/3.0 * np.exp(-(X+1)**2 - Y**2)


def pad(A, width=1, edges="all", value=0.0):
    """ Apply padding to a 2D array.
            *A*         :   array to pad
            *width*     :   thickness of padding
            *edges*     :   "all", "left", "right", "top", "bottom"
            *value"     :   0.0, value to pad with
    """
    ny = A.shape[0]
    nx = A.shape[1]

    if edges == "all":
        edges = ["left", "right", "bottom", "top"]

    if "top" == edges or "top" in edges:
        ny += width
        y0 = width
    else:
        y0 = 0
    if "bottom" == edges or "bottom" in edges:
        ny += width
        yf = -width
    else:
        yf = None
    if "left" == edges or "left" in edges:
        nx += width
        x0 = width
    else:
        x0 = 0
    if "right" == edges or "right" in edges:
        nx += width
        xf = -width
    else:
        xf = None

    B = value * np.ones([ny, nx])
    B[y0:yf, x0:xf] = A

    return B


def slope(D, res=(30.0, 30.0)):
    """ Return the scalar slope at each pixel. Use the neighbourhood
    method.
    http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=How%20Slope%20works
    """
    dx = res[0]
    dy = res[1]
    Ddx = ((2 * D[1:-1,2:] + D[:-2,2:] + D[2:,2:]) -
           (2 * D[1:-1,:-2] + D[:-2,:-2] + D[2:,:-2])) / (8.0 * dx)
    Ddy = ((2 * D[2:,1:-1] + D[2:,2:] + D[2:,:-2]) -
           (2 * D[:-2,1:-1] + D[:-2,:-2] + D[:-2,2:])) / (8.0 * dy)

    return pad(np.sqrt(Ddx*Ddx + Ddy*Ddy))


def aspect(D, res=(30.0, 30.0)):
    """ Return the slope aspect for each pixel.
    http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=How%20Aspect%20works
    """
    dx = res[0]
    dy = res[1]
    Ddx = ((2 * D[1:-1,2:] + D[:-2,2:] + D[2:,2:]) -
           (2 * D[1:-1,:-2] + D[:-2,:-2] + D[2:,:-2])) / (8.0 * dx)
    Ddy = ((2 * D[2:,1:-1] + D[2:,2:] + D[2:,:-2]) -
           (2 * D[:-2,1:-1] + D[:-2,:-2] + D[:-2,2:])) / (8.0 * dy)

    return pad(np.arctan2(Ddy, -Ddx))


def grad(D, res=(30.0, 30.0)):
    """ Computes the gradient of potential D. Return a tuple (dx, dy).
    """
    dx = res[0]
    dy = res[1]
    Ddx = ((2 * D[1:-1,2:] + D[:-2,2:] + D[2:,2:]) -
           (2 * D[1:-1,:-2] + D[:-2,:-2] + D[2:,:-2])) / (8.0 * dx)
    Ddy = ((2 * D[2:,1:-1] + D[2:,2:] + D[2:,:-2]) -
           (2 * D[:-2,1:-1] + D[:-2,:-2] + D[:-2,2:])) / (8.0 * dy)
    return pad(Ddx, width=1, edges='all'), pad(Ddy, width=1, edges='all')


def div(U, V, res=(30.0, 30.0)):
    """ Calculate the divergence of a vector field. """
    dUdx = (U[:,2:] - U[:,:-2]) / (2.0*res[0])
    dVdy = (V[2:,:] - V[:-2,:]) / (2.0*res[1])
    divergence = pad(dUdx, width=1, edges=('left', 'right')) \
               + pad(dVdy, width=1, edges=('top', 'bottom'))
    return divergence


def normed_vector_field(D):
    """ Computes a U,V vector field of potential D. Scalar components of
    U,V are normalized to max(|U, V|).
    """
    Ddx = ((2 * D[1:-1,2:] + D[:-2,2:] + D[2:,2:]) -
           (2 * D[1:-1,:-2] + D[:-2,:-2] + D[2:,:-2])) / 8.0
    Ddy = ((2 * D[2:,1:-1] + D[2:,2:] + D[2:,:-2]) -
           (2 * D[:-2,1:-1] + D[:-2,:-2] + D[:-2,2:])) / 8.0
    M = np.sqrt(Ddx**2 + Ddy**2)
    U = Ddx / M[np.isnan(M)==False].max()
    V = Ddy / M[np.isnan(M)==False].max()
    return pad(U), pad(V)


def hillshade(D, res=(30.0, 30.0), bearing=330.0, azimuth=60.0):
    """ Return a hillshade raster for field D. Testing version. """
    dx, dy = grad(D, res=res)
    u = np.array((res[0] * np.ones_like(dx), np.zeros_like(dx), dx))
    v = np.array((np.zeros_like(dy), res[1] * np.ones_like(dy), dy))
    w = np.cross(u, v, axisa=0, axisb=0)
    wunit = w / np.atleast_3d(np.sqrt(np.sum(w**2, axis=-1)))
    s = np.array((np.cos(bearing*pi/180.0),
                  np.sin(bearing*pi/180.0),
                  np.sin(azimuth*pi/180.0)))
    smat = s*np.ones([wunit.shape[0], wunit.shape[1], 3])
    dprod = (wunit*smat).sum(axis=-1)
    return dprod.T


def neighbours_of(a):
    """ For a (z,i,j) point `a`, return the neighbouring indices. """
    _, i, j = a
    return ((i-1, j-1), (i, j-1), (i+1, j-1),
            (i-1, j), (i+1, j),
            (i-1, j+1), (i, j+1), (i+1, j+1))


def fill_sinks(Z):
    """ Fill sinks in a DEM following the algorithm of Wang and Liu
    (2006).

        *Z*    :    2d array of elevation or potential data
                    (must not contain NaN!)

    Wang, L. and Liu, H. An efficient method for identifying and filling
    surface depressions in digital elevation models for hydrologic
    analysis and modelling. International Journal of Geographical
    Information Science, 20:2 (2006).
    """

    # Initialize SPILL and CLOSED
    SPILL = Z.copy()
    CLOSED = np.zeros_like(Z)
    OPEN = []

    # Get the boundary cells
    ny, nx = (Z.shape[0], Z.shape[1])
    B = [(i, j) for i in range(ny) for j in (0, nx-1)]
    B.extend([(i, j) for i in (0, ny-1) for j in range(nx)])

    # Get z along the boundary
    for b in B:
        SPILL[b] = Z[b]
        OPEN.append((Z[b], b[0], b[1]))

    while len(OPEN) > 0:

        OPEN.sort()
        c = OPEN.pop(0)
        CLOSED[c[1:]] = 1

        for n in neighbours_of(c):
            if (n[0]<ny) and (n[0]>=0) and (n[1]<nx) and (n[1]>=0):
                Zn = Z[n]
                if CLOSED[n] == 0 and (Zn, n[0], n[1]) not in OPEN:
                    SPILL[n] = max(Zn, SPILL[c[1], c[2]])
                    OPEN.append((Zn, n[0], n[1]))

    return SPILL


def viewshed(D, i, j, r=-1):
    """ Return the viewshed on *D* at (i,j).
    """
    raise NotImplementedError
