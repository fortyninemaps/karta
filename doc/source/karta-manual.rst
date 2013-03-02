.. karta documentation master file, created by
   sphinx-quickstart on Sat Feb 16 10:13:48 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to karta's documentation!
=====================================

Contents:

.. toctree::
   :maxdepth: 2

Introduction
------------

*Karta* contains a collection of loosely-related Python modules for performing
lightweight geospatial data analysis. *Karta* replicates the functionality of
GDAL, OGR, and GSTAT in some cases, but does so in a minimal package with few
dependencies, and attempts to solve problems in ways that would be considered
'pythonic'. Although these other packages may have advantages in performance or
capability, *Karta* is designed to be portable require as little extra software
infrastructure as possible.

*Karta* should be considered experimental, and no assurances are provided about
the correctness of the code or results derived using it.


Package reference
-----------------

Vector package
~~~~~~~~~~~~~~

.. automodule:: karta.vector.guppy
    :members:

.. automodule:: karta.vector.geojson
    :members:

.. automodule:: karta.vector.xyfile
    :members:

Raster package
~~~~~~~~~~~~~~

.. automodule:: karta.raster.raster
    :members:

.. automodule:: karta.raster.aaigrid
    :members:

.. automodule:: karta.raster.grid
    :members:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

