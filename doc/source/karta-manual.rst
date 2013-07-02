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

*Karta* is a Leatherman for geographic analyses in Python. *Karta* provides an
interface for solving problems in Python. To this end, it provides a simple and
clean vector and raster data types, a selection of analysis functions, the
ability to read and write a useful subset of formats, and close interoperability
with *numpy*.

Goals of *Karta* include providing a simple, lightweight, and fast set of tools
useful for everyday spatial analysis. *Karta* should be considered a work in
progress.

Formats
~~~~~~~

*Karta* attempts to provide a basic working interface to several of common file formats.
Currently partially supported are:

- vector
    - ASCII tables (XYZ) (r,w)
    - GeoJSON (r,w)
    - VTK (w)
    - ESRI Shapefiles via pyshp (r,w)
- raster
    - ESRI ASCII Grids (r,w)
    - GeoTiff (WIP)

Dependencies
~~~~~~~~~~~~

- Python 2.x
- numpy
- scipy (optional)

**CYTHON**

Cython is an optional dependency used to speed up select functions. To compile the
Cython-enabled sub-modules, run:

    setup.py build_ext --inplace

In general, enhanced-performance functions will then be called automatically when
available, otherwise *Karta* will fall back to numpy and pure-Python versions.

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

