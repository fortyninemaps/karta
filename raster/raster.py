"""
Raster functions
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


def slope(D, res=[30.0, 30.0]):
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


def aspect(D, res=[30.0, 30.0]):
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


def vector_field(D):
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


def fill_sinks(Z):
    """ Fill sinks in a DEM following the algorithm of Wang and Liu
    (2006).

    Wang, L. and Liu, H. An efficient method for identifying and filling
    surface depressions in digital elevation models for hydrologic
    analysis and modelling. International Journal of Geographical
    Information Science, 20:2 (2006).
    """

    def neighbours_of(a):
        i, j, z = a
        return ((i-1, j-1), (i, j-1), (i+1, j-1),
                (i-1, j), (i+1, j),
                (i-1, j+1), (i, j+1), (i+1, j+1))

    # Initialize SPILL and CLOSED
    SPILL = Z.copy()
    CLOSED = np.zeros_like(Z)
    OPEN = []

    # Get the boundary cells
    ny, nx = Z.shape
    B = [(i, j) for i in range(ny) for j in (0, nx-1)]
    B.extend([(i, j) for i in (0, ny-1) for j in range(nx)])

    # Get z along the boundary
    for b in B:
        SPILL[b] = Z[b]
        OPEN.append((b[0], b[1], Z[b]))

    while len(OPEN) > 0:

        OPEN.sort(key=lambda a: a[2], reverse=True)
        c = OPEN.pop()
        CLOSED[c[:2]] = 1

        N = neighbours_of(c)
        for n in N:
            if (n[0]>=ny or n[0]<0 or n[1]>=nx or n[1]<0) is False:
                Zn = Z[n]
                if CLOSED[n] == 0 and (n[0], n[1], Zn) not in OPEN:
                    SPILL[n] = max(Zn, SPILL[c[0], c[1]])
                    OPEN.append((n[0], n[1], Zn))

    return SPILL


def viewshed(D, i, j, r=-1):
    """ Return the viewshed on *D* at (i,j).
    """
