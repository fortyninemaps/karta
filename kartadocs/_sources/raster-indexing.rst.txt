Raster indexing
===============

Raster indexing is handled by a ``BandIndexer``.

Grid slicing
------------

The ``RegularGrid`` object represents a multiband raster image. The underlying
data can be extracted as a numpy array using slicing syntax. In general, a slice
preserves a dimension, while an integer collapses it. To illustrate:

::

    grid[:,:,:]     # Returns a three-dimensional numpy array (nrow, ncol, nband)
    grid[:,:]       # Shorthand for the above

    grid[:,:,0]     # Returns a two-dimensional numpy array representing the
                    # first band (nrow, ncol)

    grid[:,:,:2]    # Returns a three dimensional numpy array with the first two
                    # bands (nrow, ncol, 2)

    grid[:,10]      # Returns all band values for the tenth column (nrow, nband)
    grid[:,10,0]    # Returns the first band values for the tenth column (nrow,)

Grid masking
------------

It is possible to extract grid values using numpy arrays as masks. This works
similarly to indexing a numpy array with shape (nrow, ncol, nband).

::

    mask.shape      # (nrow, ncol, nband)
    sum(mask)       # n

    grid[mask3]     # an (n,) array

It's also possible to index with a two-dimensional array. Unlike in numpy, the
mask is broadcast along the band dimension.

::

    mask.shape      # (nrow, ncol)
    sum(mask)       # n

    grid[mask2]     # an (n, nband) array

Warning: as currently implemented, masking with an array is syntactic sugar
around unpacking the entire array and indexing with numpy array operations. As a
result, array masking generally results in the entire grid being decompressed or
read into memory, depending on the ``Band`` backend.

