Karta - simple geospatial analysis in Python
--------------------------------------------

*Karta* is a Leatherman for geographic analyses in Python. *Karta* provides an API for
solving problems in Python. To this end, it provides a simple and clean interface to both
vector and raster data types, the ability to read and write a useful subset of formats,
and close interoperability with *numpy*.

Goals of *Karta* include providing a simple, lightweight, and fast set of tools useful for
everyday spatial analysis. *Karta* should be considered a work in progress.

**Curently working on:**
- projection handling
- refactoring raster class hierarchy
- native GeoTiff support

**Thinking about:**
- memory efficiency

##CONTENTS

- vector
    - guppy :       Basic vector geometry classes (e.g. `Point`, `Multipoint`, `Line`, `Polygon`)
    - gpx_parser :  Parser for GPX files exported from GPS devices
    - geojson :     Classes and functions for handling GeoJSON files
    - vtk :         XML-based VTK interface
    - shp_funcs :   Shapefile-to-guppy conversions
    - stats :       Geostatistical functions
    - xyfile :      ASCII table functions

- raster
    - grid :        Basic Grid types, including `StructuredGrid` and `RegularGrid`
    - aaigrid :     Grid subclass specifically for reading, writing, and manipulating ESRI ASCII grids
    - flow :        Stream flow functions
    - raster :      General purpose raster functions
    - streamline :  Streamline calculation

- tests : unit tests


##FORMATS

*Karta* attempts to provide a basic working interface to several of common file formats.
Currently partially supported are:

- vector
    - ASCII tables (XYZ) (r,w)
    - GeoJSON (r,w)
    - VTK (w)
    - ESRI Shapefiles (r,w)
- raster
    - ESRI ASCII Grids (r,w)
    - GeoTiff (WIP)

##CYTHON

Cython is an optional dependency than can make select functions faster. To
compile the Cython-enabled sub-modules, run:

    setup.py build_ext

In general, enhanced-performance functions will then be called automatically
where available.



##LICENSE

This software is provided under the MIT license.

The vector module contains a snapshot of the pyshp shapefile module, available
at (http://code.google.com/p/pyshp/). This module is also available under the
MIT license.

###MIT License:

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

