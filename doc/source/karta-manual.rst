.. karta documentation master file, created by
   sphinx-quickstart on Sat Feb 16 10:13:48 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Karta
=====

Contents:

.. toctree::
   :maxdepth: 2

   installation
   tutorial
   reference

*Karta* provides a simple and fast framework for spatial analysis in Python.

The package provides clean vector and raster data types that are geographical
coordinate system-aware, a selection of geographical analysis methods, and the
ability to read and write several formats, including GeoJSON, shapefiles, and
ESRI ASCII.

*Karta* works with Python 2.6-2.7 and Python 3.3+. Suggestions, bug reports,
testcases, and pull requests are welcome, particularly to improve file format
support and test coverage.

See the tutorial_, or read the `API documentation`_.

.. _tutorial: tutorial.html
.. _API documentation: reference.html

Data Formats
------------

*Karta* provides a basic working interface to several of common file
formats. Currently supported are:

-  vector

   -  GeoJSON (r,w)
   -  ESRI Shapefiles via pyshp (r,w)
   -  ASCII tables (XYZ) (r,w)
   -  GPS eXchange (GPX) (r,w)

-  raster

   -  ESRI ASCII Grid (r,w)
   -  GeoTiff (through GDAL) (r)
   -  USGS DEM (WIP)

Examples
--------

`Sharing data via __geo_interface__ <geointerface.html>`_

*(work in progress)*


License
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

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

