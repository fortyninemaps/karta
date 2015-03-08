# Karta

[![Build Status](https://travis-ci.org/njwilson23/karta.svg?branch=master)](https://travis-ci.org/njwilson23/karta)

*Karta* is a simple and fast framework for spatial analysis in Python.

It contains clean vector and raster data types that are coordinate system-aware,
implementations of frequently-used of geospatial analysis methods, and the
read/write interfaces to several formats, including GeoJSON, shapefiles, and
ESRI ASCII.

*Karta* works with Python 2.6-2.7 and Python 3.3+. Suggestions, bug reports,
test cases, and pull requests are welcome.

## DOCUMENTATION

See the [online manual](http://www.ironicmtn.com/kartadocs/karta-manual.html),
the [tutorial](http://www.ironicmtn.com/kartadocs/tutorial.html), or read the
[API documentation](http://www.ironicmtn.com/kartadocs/reference.html).

The manual can also be built offline using Sphinx by running `make` from the
`doc/` subdirectory. The documentation is built from source code docstrings and
the example IPython notebooks, which are also reproduced in the
[Wiki](https://github.com/njwilson23/karta/wiki/Tutorial). Building the
documentation requires [Sphinx](http://sphinx-doc.org/),
[alabaster](https://github.com/bitprophet/alabaster) and
[numpydoc](https://github.com/numpy/numpydoc).

## PACKAGE OVERVIEW

- **karta.crs**: framework for coordinate reference systems and geodetic
  calculations

- **karta.vector.geometry**: geometry classes `Point`, `Multipoint`, `Line`, and
  `Polygon` with associated methods such as length, area, intersections,
  membership testing, convex hulls, and affine transformations

- **karta.raster.grid**: `Grid` classes including `RegularGrid` class
  (supporting CRS-aware clipping, sampling, profiling along vector tracks), and
  experimental `WarpedGrid`

- **karta.tests**: unit tests, to be run with `python tests/runtests.py`

## FORMATS

*Karta* provides a basic working interface to several of common file formats.
Currently partially-supported are:

- vector
    - GeoJSON (r,w)
    - ESRI Shapefiles (via pyshp) (r,w)
    - GPS eXchange (GPX) (r,w)
    - ASCII tables (XYZ) (r,w)
- raster
    - ESRI ASCII Grid (r,w)
    - GeoTiff (via GDAL) (r,w)
    - USGS DEM (WIP)

*Karta* implements [Python
`__geo_interface__`](https://gist.github.com/sgillies/2217756) for vector
geometries. External modules that also implement `__geo_interface__` (e.g.
[shapely](https://github.com/Toblerity/Shapely),
[fastkml](https://fastkml.readthedocs.org/en/latest/)) can be directly used to
read and write additional formats.

## INSTALLATION

The easiest way to install in production is to use `pip`. Installation requires
a version of `setuptools>=0.7.0`.

    pip install -U setuptools

To install the latest release from PyPI, run

    pip install karta

To build from source,

    git clone https://github.com/njwilson23/karta.git
    pip install -r karta/requirements.txt
    pip install karta/

## DEPENDENCIES

### Required

- Python 2.6+ or Python 3.3+
- numpy
- pyshp
- pyproj

### Recommended

- cython
- gdal (for geotiff I/O)

### Optional

- scipy

When installing from PyPI, Cython-compiled C source code is provided and will be
automatically compiled to improve performance if a suitable C compiler is
available.

## LICENSE

This software is provided under the MIT license.

### MIT License:

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

