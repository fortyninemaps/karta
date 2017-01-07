"""
2D raster functions
"""

import numpy as np
from .grid import RegularGrid

def _slope(D, res=(1.0, 1.0)):
    """ Return the scalar slope at each pixel using the neighbourhood method.
    http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=How%20Slope%20works
    """
    dx, dy = res
    Ddx = ((2 * D[1:-1,2:] + D[:-2,2:] + D[2:,2:]) -
           (2 * D[1:-1,:-2] + D[:-2,:-2] + D[2:,:-2])) / (8.0 * dx)
    Ddy = ((2 * D[2:,1:-1] + D[2:,2:] + D[2:,:-2]) -
           (2 * D[:-2,1:-1] + D[:-2,:-2] + D[:-2,2:])) / (8.0 * dy)
    return np.pad(np.sqrt(Ddx*Ddx + Ddy*Ddy), ((1, 1), (1, 1)), "reflect",
                  reflect_type="odd")

def slope(grid, band=0):
    """ Return the scalar slope at each pixel using the neighbourhood method.

    Parameters
    ----------
    grid: RegularGrid
    band: int, optional
        band to compute hillshade for (default 0)

    Returns
    -------
    RegularGrid

    Notes
    -----
    http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=How%20Slope%20works
    """
    if grid.skew != (0, 0):
        raise NotImplementedError("slope calculations not implemented on skewed grids")
    D = np.where(grid.data_mask, grid[:,:,band].astype(np.float32), np.nan)
    return RegularGrid(grid.transform, _slope(D, grid.resolution),
                       crs=grid.crs, nodata_value=np.nan)

def _aspect(D, res=(1.0, 1.0)):
    """ Return the slope aspect for each pixel.
    http://webhelp.esri.com/arcgisdesktop/9.2/index.cfm?TopicName=How%20Aspect%20works
    """
    Ddx = ((2 * D[1:-1,2:] + D[:-2,2:] + D[2:,2:]) -
           (2 * D[1:-1,:-2] + D[:-2,:-2] + D[2:,:-2])) / (8.0 * res[0])
    Ddy = ((2 * D[2:,1:-1] + D[2:,2:] + D[2:,:-2]) -
           (2 * D[:-2,1:-1] + D[:-2,:-2] + D[:-2,2:])) / (8.0 * res[1])
    return np.pad(np.arctan2(Ddy, -Ddx), ((1, 1), (1, 1)), "constant",
                  constant_values=(np.nan,))

def aspect(grid, band=0):
    """ Compute grid aspect.

    Parameters
    ----------
    grid: RegularGrid
    band: int, optional
        band to compute hillshade for (default 0)

    Returns
    -------
    RegularGrid
    """
    if grid.skew != (0, 0):
        raise NotImplementedError("aspect calculations not implemented on skewed grids")
    D = np.where(grid.data_mask, grid[:,:,band].astype(np.float32), np.nan)
    return RegularGrid(grid.transform, _aspect(D, grid.resolution),
                       crs=grid.crs, nodata_value=np.nan)

def _grad(D, res=(1.0, 1.0)):
    """ Computes the gradient of potential D. Return a tuple (dx, dy).
    """
    Ddx = ((2 * D[1:-1,2:] + D[:-2,2:] + D[2:,2:]) -
           (2 * D[1:-1,:-2] + D[:-2,:-2] + D[2:,:-2])) / (8.0 * res[0])
    Ddy = ((2 * D[2:,1:-1] + D[2:,2:] + D[2:,:-2]) -
           (2 * D[:-2,1:-1] + D[:-2,:-2] + D[:-2,2:])) / (8.0 * res[1])
    return (np.pad(Ddx, ((1, 1), (1, 1)), "constant", constant_values=(np.nan,)),
            np.pad(Ddy, ((1, 1), (1, 1)), "constant", constant_values=(np.nan,)))

def gradient(grid, band=0):
    """ Compute gradient field from a grid.

    Parameters
    ----------
    grid: RegularGrid
    band: int, optional (default 0)

    Returns
    -------
    (RegularGrid, RegularGrid)
    """
    if grid.skew != (0, 0):
        raise NotImplementedError("gradient calculations not implemented on skewed grids")
    D = np.where(grid.data_mask, grid[:,:,band].astype(np.float32), np.nan)
    dDdx, dDdy = _grad(D, grid.resolution)
    return (RegularGrid(grid.transform, dDdx, crs=grid.crs, nodata_value=np.nan),
            RegularGrid(grid.transform, dDdy, crs=grid.crs, nodata_value=np.nan))

def _div(U, V, res=(1.0, 1.0)):
    """ Calculate the divergence of a vector field. """
    dUdx = (U[:,2:] - U[:,:-2]) / (2.0*res[0])
    dVdy = (V[2:,:] - V[:-2,:]) / (2.0*res[1])
    divergence = np.pad(dUdx, ((0, 0), (1, 1)), "constant", constant_values=(np.nan,)) \
               + np.pad(dVdy, ((1, 1), (0, 0)), "constant", constant_values=(np.nan,))
    return divergence

def divergence(grid, bands=(0, 1)):
    """ Compute divergence from a grid.

    Parameters
    ----------
    grid: RegularGrid
    band: tuple of ints, optional (default (0, 1))

    Returns
    -------
    RegularGrid
    """
    if grid.skew != (0, 0):
        raise NotImplementedError("divergence calculations not implemented on skewed grids")
    return RegularGrid(grid.transform,
                       _div(grid[:,:,bands[0]],
                            grid[:,:,bands[1]],
                            grid.resolution),
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
    return (np.pad(U, ((1, 1), (1, 1)), "constant", constant_values=(np.nan,)),
            np.pad(V, ((1, 1), (1, 1)), "constant", constant_values=(np.nan,)))

def normed_potential_vectors(grid, band=0):
    """ Computes a U,V vector field of a potential grid. Scalar components of
    U,V are normalized to max(|U, V|).

    Parameters
    ----------
    grid: RegularGrid
    band: int, optional (default 0)

    Returns
    -------
    (RegularGrid, RegularGrid)
    """
    if grid.skew != (0, 0):
        raise NotImplementedError("potential vector calculations not implemented on skewed grids")
    D = np.where(grid.data_mask, grid[:,:,band].astype(np.float32), np.nan)
    u, v = _normed_potential_vectors(D, res=grid.resolution)
    return (RegularGrid(grid.transform, u, crs=grid.crs, nodata_value=np.nan),
            RegularGrid(grid.transform, v, crs=grid.crs, nodata_value=np.nan))

def hillshade(grid, azimuth=330.0, elevation=60.0, band=0):
    """ Return a hill-shaded version of *grid*.

    Parameters
    ----------
    grid: RegularGrid
    azimuth: float, optional
        direction of light source (default 330.0)
    elevation: float, optional
        height of light source (default 60.0)
    band: int, optional
        band to compute hillshade for (default 0)

    Returns
    -------
    RegularGrid

    Notes
    -----
    Currently assumes orthogonal coordinates.
    """
    dxgrid, dygrid = gradient(grid, band=band)
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

