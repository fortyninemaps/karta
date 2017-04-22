Raster indexing
===============

Raster indexing is handled by a ``BandIndexer``.

Grid slicing
------------

The ``RegularGrid`` object represent a multiband raster image. The underlying data
can be extracted as a numpy array using slicing syntax. In general, a slice
preserves a dimension, while an integer collapses it. To illustrate:

::

    grid[:,:,:]     # Returns a three-dimensional numpy array (ny, nx, nbands)
    grid[:,:]       # Shorthand for the above

    grid[:,:,0]     # Returns a two-dimensional numpy array representing the
                    # first band (ny, nx)

    grid[:,:,:2]    # Returns a three dimensional numpy array with the first two
                    # bands (ny, nx, 2)

    grid[:,10]      # Returns all band values for the tenth column (ny, nbands)
    grid[:,10,0]    # Returns the first band values for the tenth column (ny,)

Grid masking
------------

It is possible to extract grid values using numpy arrays as masks. Unlike numpy
arrays, the result has an extra dimension storing separate bands.

::

    # `mask` is a two dimensional array with *n* True values

    grid[mask]      # Returns an (n, nbands) array

    grid[mask,0]    # Returns an (n,) array representing the first band

Finally, a three-dimensional mask can be used, which collapses bands

::

    # `mask` is a three dimensional array with *n* True values

    grid[mask3]     # Returns an (n,) array

    grid[mask3,0]   # Raises an IndexError


Warning: as currently implemented, masking with an array may result in the
entire grid being decompressed or read into memory, depending on the ``Band``
backend.

