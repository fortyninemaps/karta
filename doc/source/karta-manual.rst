.. geo_tools documentation master file, created by
   sphinx-quickstart on Sat Feb 16 10:13:48 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to geo_tools's documentation!
=====================================

Contents:

.. toctree::
   :maxdepth: 2

Introduction
------------

*geo_tools* contains a collection of loosely-related Python modules for
performing lightweight geospatial data analysis. *geo_tools* replicates the
functionality of GDAL, OGR, and GSTAT in some cases, but does so in a minimal
package with few dependencies, and attempts to solve problems in ways that
would be considered 'pythonic'. Although these other packages may have
advantages in performance or capability, *geo_tools* is designed to be portable
require as little extra software infrastructure as possible.

*geo_tools* should be considered experimental, and no assurances are provided
about the correctness of the code or results derived using it.


Package reference
-----------------

Vector package
~~~~~~~~~~~~~~

.. automodule:: geo_tools.vector.guppy
    :members:

.. automodule:: geo_tools.vector.json
    :members:

.. automodule:: geo_tools.vector.xy
    :members:

Raster package
~~~~~~~~~~~~~~

.. automodule:: geo_tools.raster.raster
    :members:

.. automodule:: geo_tools.raster.aaigrid
    :members:

.. automodule:: geo_tools.raster.grid
    :members:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

