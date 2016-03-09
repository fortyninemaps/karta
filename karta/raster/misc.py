"""
2D raster functions
"""

import numpy as np
from .grid import RegularGrid

def witch_of_agnesi(nx=100, ny=100, a=4.0):
    """ Return a raster field defined by the equation
                    8 a^3
            Z = --------------
                 d^2 + 4 a*2
    where d is the distance from the center.
    """
    xc = int(np.floor(nx / 2.0))
    yc = int(np.floor(ny / 2.0))
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

def _slope(D, res=(1.0, 1.0)):
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

    return pad(np.sqrt(Ddx*Ddx + Ddy*Ddy), value=np.nan)

def slope(grid):
    """ Return the scalar slope at each pixel. Use the neighbourhood
    method.
    http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=How%20Slope%20works
    """
    if grid.transform[4:] != (0, 0):
        raise NotImplementedError("slope calculations no implemented on skewed grids")
    D = np.where(grid.data_mask, grid.values.astype(np.float32), np.nan)
    return RegularGrid(grid.transform, _slope(D, grid.resolution),
                       crs=grid.crs, nodata_value=np.nan)

def _aspect(D, res=(1.0, 1.0)):
    """ Return the slope aspect for each pixel.
    http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=How%20Aspect%20works
    """
    dx = res[0]
    dy = res[1]
    Ddx = ((2 * D[1:-1,2:] + D[:-2,2:] + D[2:,2:]) -
           (2 * D[1:-1,:-2] + D[:-2,:-2] + D[2:,:-2])) / (8.0 * dx)
    Ddy = ((2 * D[2:,1:-1] + D[2:,2:] + D[2:,:-2]) -
           (2 * D[:-2,1:-1] + D[:-2,:-2] + D[:-2,2:])) / (8.0 * dy)

    return pad(np.arctan2(Ddy, -Ddx), value=np.nan)

def aspect(grid):
    if grid.transform[4:] != (0, 0):
        raise NotImplementedError("aspect calculations no implemented on skewed grids")
    D = np.where(grid.data_mask, grid.values.astype(np.float32), np.nan)
    return RegularGrid(grid.transform, _aspect(D, grid.resolution),
                       crs=grid.crs, nodata_value=np.nan)

def _grad(D, res=(1.0, 1.0)):
    """ Computes the gradient of potential D. Return a tuple (dx, dy).
    """
    dx = res[0]
    dy = res[1]
    Ddx = ((2 * D[1:-1,2:] + D[:-2,2:] + D[2:,2:]) -
           (2 * D[1:-1,:-2] + D[:-2,:-2] + D[2:,:-2])) / (8.0 * dx)
    Ddy = ((2 * D[2:,1:-1] + D[2:,2:] + D[2:,:-2]) -
           (2 * D[:-2,1:-1] + D[:-2,:-2] + D[:-2,2:])) / (8.0 * dy)
    return (pad(Ddx, width=1, edges='all', value=np.nan),
            pad(Ddy, width=1, edges='all', value=np.nan))

def gradient(grid):
    if grid.transform[4:] != (0, 0):
        raise NotImplementedError("gradient calculations no implemented on skewed grids")
    D = np.where(grid.data_mask, grid.values.astype(np.float32), np.nan)
    dDdx, dDdy = _grad(D, grid.resolution)
    return (RegularGrid(grid.transform, dDdx, crs=grid.crs, nodata_value=np.nan),
            RegularGrid(grid.transform, dDdy, crs=grid.crs, nodata_value=np.nan))

def _div(U, V, res=(1.0, 1.0)):
    """ Calculate the divergence of a vector field. """
    dUdx = (U[:,2:] - U[:,:-2]) / (2.0*res[0])
    dVdy = (V[2:,:] - V[:-2,:]) / (2.0*res[1])
    divergence = pad(dUdx, width=1, edges=('left', 'right'), value=np.nan) \
               + pad(dVdy, width=1, edges=('top', 'bottom'), value=np.nan)
    return divergence

def divergence(grid):
    if grid.transform[4:] != (0, 0):
        raise NotImplementedError("divergence calculations no implemented on skewed grids")
    D = np.where(grid.data_mask, grid.values.astype(np.float32), np.nan)
    return RegularGrid(grid.transform, _div(D, grid.resolution),
                       crs=grid.crs, nodata_value=np.nan)

def _normed_potential_vectors(D, res=(1.0, 1.0)):
    """ Computes a U,V vector field of potential D. Scalar components of
    U,V are normalized to max(|U, V|).
    """
    dx, dy = res
    Ddx = ((2 * D[1:-1,2:] + D[:-2,2:] + D[2:,2:]) -
           (2 * D[1:-1,:-2] + D[:-2,:-2] + D[2:,:-2])) / (8.0*dx)
    Ddy = ((2 * D[2:,1:-1] + D[2:,2:] + D[2:,:-2]) -
           (2 * D[:-2,1:-1] + D[:-2,:-2] + D[:-2,2:])) / (8.0*dy)
    M = np.sqrt(Ddx**2 + Ddy**2)
    U = Ddx / M[np.isnan(M)==False].max()
    V = Ddy / M[np.isnan(M)==False].max()
    return pad(U, value=np.nan), pad(V, value=np.nan)

def normed_potential_vectors(grid):
    if grid.transform[4:] != (0, 0):
        raise NotImplementedError("vector field calculations no implemented on skewed grids")
    D = np.where(grid.data_mask, grid.values.astype(np.float32), np.nan)
    u, v = _normed_potential_vectors(D, res=grid.resolution)
    return (RegularGrid(grid.transform, u, crs=grid.crs, nodata_value=np.nan),
            RegularGrid(grid.transform, v, crs=grid.crs, nodata_value=np.nan))

def _hillshade(D, res=(1.0, 1.0), bearing=330.0, azimuth=60.0):
    """ Return a hillshade raster for field D. Testing version. """
    dx, dy = grad(D, res=res)
    u = np.array((res[0] * np.ones_like(dx), np.zeros_like(dx), dx))
    v = np.array((np.zeros_like(dy), res[1] * np.ones_like(dy), dy))
    w = np.cross(u, v, axisa=0, axisb=0)
    wunit = w / np.atleast_3d(np.sqrt(np.sum(w**2, axis=-1)))
    pi = np.pi
    s = np.array((np.cos(bearing*pi/180.0),
                  np.sin(bearing*pi/180.0),
                  np.sin(azimuth*pi/180.0)))
    smat = s*np.ones([wunit.shape[0], wunit.shape[1], 3])
    dprod = (wunit*smat).sum(axis=-1)
    return dprod

def hillshade(grid, **kw):
    """ Return a hill-shaded version of *grid*.

    Arguments
    ---------
    grid: RegularGrid instance
    bearing: float, bearing of light source (default 330.0)
    azimuth: float, azimuth of light source (default 60.0)

    Note: Currently assumes orthogonal coordinates.
    """
    hs = RegularGrid(grid.transform,
                     values=_hillshade(grid.values, res=grid.resolution, **kw),
                     crs=grid.crs)
    q = np.percentile(hs.values[~np.isnan(hs.values)], [2, 98])
    np.clip(hs.values, q[0], q[1], out=hs.values)
    return hs

def viewshed(D, i, j, r=-1):
    """ Return the viewshed on *D* at (i,j).
    """
    raise NotImplementedError()
