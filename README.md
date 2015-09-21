# Karta

[![Build Status](https://travis-ci.org/njwilson23/karta.svg?branch=master)](https://travis-ci.org/njwilson23/karta)

*Karta* is a simple and fast framework for spatial analysis in Python.

Create vector geometries:

```python
point = Point((-130.0, 52.0), crs=LonLatWGS84)

line = read_geojson("linedata.json")

polygon = Polygon([(-515005.78, -1301130.53),
                   (-579174.89, -1282271.94),
                   (-542977.83, -1221147.82),
                   (-437864.05, -1251641.55),
                   (-438160.72, -1252421.48),
                   (-437961.28, -1285314.00)],
                   crs=NSIDCNorth)
```
Perform simple queries:
```python
point2 = Point((-25.0, 48.0), crs=LonLatWGS84)
point.distance(point2)          # Distance in geographical units

line.intersects(polygon)        # True or False

ch = polygon.convex_hull()      # Returns a new polygon
ch.to_shapefile("poly.shp")
```
Work with raster data:
```python
grid = read_gtiff("landsat_scene.tif")  # Leverages GDAL

grid.profile(line)              # Collect data along a line

grid.resample(500.0, 500.0)     # Return a grid resampled at a new resolution
```

The latest release is on PyPI (see [Installation](#installation)).
*Karta* works with Python 2.7 and Python 3.3+. Suggestions, bug reports, test
cases, and pull requests are welcome.

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
Currently implemented are:

- vector
    - GeoJSON (r,w)
    - ESRI Shapefiles (via pyshp) (r,w)
    - GPS eXchange (GPX) (r,w)
- raster
    - ESRI ASCII Grid (r,w)
    - GeoTiff (requires GDAL) (r,w)

*Karta* implements the Python [`__geo_interface__`
attribute](https://gist.github.com/sgillies/2217756) for vector geometries. This
permits data to be exchanged between *Karta* and external modules that also
implement `__geo_interface__` (e.g.
[shapely](https://github.com/Toblerity/Shapely),
[fastkml](https://fastkml.readthedocs.org/en/latest/)).

## INSTALLATION

The easiest way to install in production is to use `pip`. Installation requires
a version of `setuptools>=17.0`.

    pip install -U setuptools

To install the latest release from PyPI, run

    pip install karta

### Building from source

Building from source requires Cython to be available.

    pip install Cython

Then, clone the repository and install,

    git clone https://github.com/njwilson23/karta.git karta
    pip install -r karta/requirements.txt
    pip install karta/

## DEPENDENCIES

### Required

- Python 2.6+ or Python 3.3+
- numpy
- pyshp
- pyproj
- C-compiler

### Recommended

- osgeom.gdal (for geotiff I/O)
- osgeo.osr (for coordinate system interchange)
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
