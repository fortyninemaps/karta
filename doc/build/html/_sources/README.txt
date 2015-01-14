Karta - tidy Python package for geospatial computation
======================================================

|Build Status|

*Karta* is a Python/Python3 package for geospatial data structures.
*Karta* serves as a Leatherman for geographic analyses. It provides
simple and clean vector and raster data types, a selection of
geographical analysis methods, and the ability to read and write several
formats, including GeoJSON, shapefiles, and ESRI ASCII.

The goals of *Karta* is to expose a simple and fast framework for
spatial analysis.

Suggestions, bug reports, testcases, and pull requests are welcome,
particularly to improve format support and test coverage.

**Current development focus:**

-  documentation
-  reprojections

DOCUMENTATION
-------------

API documentation can be built by running ``make`` from the ``doc/``
directory. The beginnings of a tutorial to introduce *Karta* are in the
`Wiki <https://github.com/njwilson23/karta/wiki/Tutorial>`__.

CONTENTS AT A GLANCE
--------------------

-  vector

   -  guppy: Vector geometry classes (e.g. ``Point``, ``Multipoint``,
      ``Line``, ``Polygon``) supporting the `Python
      \_\_geo\_interface\_\_ <https://gist.github.com/sgillies/2217756>`__
   -  gpx: GPX class for parsing and constructing GPX eXchange files
   -  geojson: Classes and functions for reading and writing GeoJSON
   -  shp\_funcs: Shapefile-to-guppy conversions through *pyshp*
      interface
   -  xyfile: ASCII table functions
   -  quadtree: QuadTree implementation
   -  stats: Geostatistical functions (experimental)

-  raster

   -  grid: Basic Grid types, including ``RegularGrid`` and
      ``WarpedGrid``
   -  raster: General purpose raster functions
   -  aaigrid: Grid subclass specifically for reading and writing ESRI
      ASCII grids
   -  streamline: Streamline calculation
   -  flow: Stream flow functions (experimental)

-  crs : ``CRS`` and ``CRSRegister`` classes used throughout *Karta*

-  tests : unit tests

FORMATS
-------

*Karta* provides a basic working interface to several of common file
formats. Currently partially-supported are:

-  vector

   -  GeoJSON (r,w)
   -  ESRI Shapefiles via pyshp (r,w)
   -  GPS eXchange (GPX) (r,w)
   -  ASCII tables (XYZ) (r,w)

-  raster

   -  ESRI ASCII Grid (r,w)
   -  GeoTiff through GDAL (r)
   -  USGS DEM (WIP)

INSTALLATION
------------

The easiest way to install is to use ``pip``. Installation requires a
version of ``setuptools>=0.7.0``.

::

    pip install -U setuptools

To install the latest release from PyPI, run

::

    pip install karta

To build from source,

::

    git clone https://github.com/njwilson23/karta.git
    pip install -r karta/requirements.txt
    pip install karta/

DEPENDENCIES
------------

Required
~~~~~~~~

-  Python 2.7+ or Python 3.3+
-  numpy
-  pyshp
-  pyproj (for geodetic calculations)

Recommended
~~~~~~~~~~~

-  cython
-  gdal (for geotiff I/O)
-  scipy

When installing from PyPI, Cython-compiled C source code is provided and
will be automatically compiled to improve performance if a suitable C
compiler is available.

TESTING
-------

To run all unit tests, execute

::

    python tests/runtests.py

LICENSE
-------

This software is provided under the MIT license.

MIT License:
~~~~~~~~~~~~

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

.. |Build Status| image:: https://travis-ci.org/njwilson23/karta.svg?branch=master
   :target: https://travis-ci.org/njwilson23/karta
