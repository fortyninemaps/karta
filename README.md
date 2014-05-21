#Karta - tidy Python package for geospatial computation

[![Build Status](https://travis-ci.org/njwilson23/karta.svg?branch=master)](https://travis-ci.org/njwilson23/karta)

*Karta* is a Leatherman for geographic analyses. *Karta* provides an interface
for solving problems in Python (2 or 3) that works nicely with existing
packages. To this end, it provides simple and clean vector and raster data
types, a selection of analysis functions, the ability to read and write a small
number of useful formats, and interoperability with *numpy*.

Goals of *Karta* include providing a simple, lightweight, and fast set of tools
useful for "everyday" spatial analysis, as well as a flexible set of
abstractions upon which to build more advanced routines. *Karta* should be
considered a work in progress.

**Current projects:**
- test coverage
- minimal system for coordinate system metadata

**Future goals:**
- native GeoTiff support?

##CONTENTS AT A GLANCE

- vector
    - guppy:        Vector geometry classes (e.g. `Point`, `Multipoint`, `Line`, `Polygon`)
    - gpx:          GPX class for parsing and constructing GPX eXchange files
    - geojson:      Classes and functions for reading and writing GeoJSON
    - vtk:          XML-based VTK interface
    - shp\_funcs:   Shapefile-to-guppy conversions through _pyshp_ interface
    - stats:        Geostatistical functions
    - xyfile:       ASCII table functions

- raster
    - grid:         Basic Grid types, including `StructuredGrid` and `RegularGrid`
    - aaigrid:      Grid subclass specifically for reading, writing, and manipulating ESRI ASCII grids
    - flow:         Stream flow functions
    - raster:       General purpose raster functions
    - streamline:   Streamline calculation

- tests : unit tests


##FORMATS

*Karta* provides a basic working interface to several of common file formats.
Currently partially-supported are:

- vector
    - ASCII tables (XYZ) (r,w)
    - GeoJSON (r,w)
    - GPS eXchange (GPX) (r,w)
    - VTK (w)
    - ESRI Shapefiles via pyshp (r,w)
- raster
    - ESRI ASCII Grid (r,w)
    - USGS DEM (WIP)

## INSTALLATION

The easiest way to install is to use `pip`.

    cd karta/
    pip install .

## DEPENDENCIES

### Required

- pyshp
- Python 2.6+ or Python 3.2+
- numpy

### Optional

- pyproj (required for geodetic calculations)
- scipy
- Cython

Cython is an optional dependency used to speed up select functions. In general,
enhanced-performance functions will then be called automatically when available,
otherwise *Karta* will fall back to numpy and pure-Python versions.

## TESTING

To run all unit tests, execute

    python tests/runtests.py

##LICENSE

This software is provided under the MIT license.

###MIT License:

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

