"""
2D raster functions
"""

import numpy as np
from .grid import RegularGrid

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

def pad(A, width=1, edges="all", value=0.0):
    """ Apply padding to a 2D array.

    Parameters
    ----------
    A : ndarray
        array to pad
    width : int
        thickness of padding
    edges : str
        one of "all", "left", "right", "top", "bottom"
    value : number
        fill value for padding

    Returns
    -------
    ndarray
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

    B = np.full((ny, nx), value, dtype=type(value))
    B[y0:yf, x0:xf] = A
    return B

def _slope(D, res=(1.0, 1.0)):
    """ Return the scalar slope at each pixel using the neighbourhood method.
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
    """ Return the scalar slope at each pixel using the neighbourhood method.
    http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=How%20Slope%20works
    """
    if grid.skew != (0, 0):
        raise NotImplementedError("slope calculations not implemented on skewed grids")
    if grid.nbands != 1:
        raise ValueError("input grid must be single-banded")
    D = np.where(grid.data_mask, grid[:,:].astype(np.float32), np.nan)
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
    if grid.skew != (0, 0):
        raise NotImplementedError("aspect calculations not implemented on skewed grids")
    if grid.nbands != 1:
        raise ValueError("input grid must be single-banded")
    D = np.where(grid.data_mask, grid[:,:].astype(np.float32), np.nan)
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
    if grid.skew != (0, 0):
        raise NotImplementedError("gradient calculations not implemented on skewed grids")
    if grid.nbands != 1:
        raise ValueError("input grid must be single-banded")
    D = np.where(grid.data_mask, grid[:,:].astype(np.float32), np.nan)
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
    if grid.skew != (0, 0):
        raise NotImplementedError("divergence calculations not implemented on skewed grids")
    if grid.nbands != 1:
        raise ValueError("input grid must be single-banded")
    D = np.where(grid.data_mask, grid[:,:].astype(np.float32), np.nan)
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
    if grid.skew != (0, 0):
        raise NotImplementedError("potential vector calculations not implemented on skewed grids")
    if grid.nbands != 1:
        raise ValueError("input grid must be single-banded")
    D = np.where(grid.data_mask, grid[:,:].astype(np.float32), np.nan)
    u, v = _normed_potential_vectors(D, res=grid.resolution)
    return (RegularGrid(grid.transform, u, crs=grid.crs, nodata_value=np.nan),
            RegularGrid(grid.transform, v, crs=grid.crs, nodata_value=np.nan))

def hillshade(grid, azimuth=330.0, elevation=60.0):
    """ Return a hill-shaded version of *grid*.

    Parameters
    ----------
    grid: RegularGrid instance
    azimuth: float, optional
        direction of light source (default 330.0)
    elevation : float, optional
        height of light source (default 60.0)

    Notes
    -----
    Currently assumes orthogonal coordinates.
    """
    if grid.nbands != 1:
        raise ValueError("input grid must be single-banded")
    dxgrid, dygrid = gradient(grid)
    dx = dxgrid[:,:]
    dy = dygrid[:,:]
    res = grid.resolution
    u = np.array((np.full_like(dx, res[0]), np.zeros_like(dx), dx))
    v = np.array((np.zeros_like(dy), np.full_like(dy, res[1]), dy))
    w = np.cross(u, v, axisa=0, axisb=0)
    wunit = w / np.atleast_3d(np.sqrt(np.sum(w**2, axis=-1)))

    s = np.array((np.cos(azimuth*np.pi/180.0),
                  np.sin(azimuth*np.pi/180.0),
                  np.sin(elevation*np.pi/180.0)))
    smat = np.full((wunit.shape[0], wunit.shape[1], 3), s, dtype=np.float64)
    dprod = (wunit*smat).sum(axis=-1)

    q = np.percentile(dprod[~np.isnan(dprod)], [2, 98])
    np.clip(dprod, q[0], q[1], out=dprod)

    return RegularGrid(grid.transform, values=dprod, crs=grid.crs)

