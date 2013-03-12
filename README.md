Karta - simple geospatial analysis in Python
--------------------------------------------

*Karta* contains a collection of loosely-related Python modules for performing
geospatial data analysis. *Karta* attempts to solve problems in ways that might
be considered 'pythonic'. To this end, it provides a simple and clean API and
close interoperability with numpy. The *Karta* package is a simple ecosystem of
modules for dealing with vector and raster data.

Goals of *Karta* include providing a simple, lightweight, and fast set of
spatially-inclined tools. *Karta* should be considered a work in progress.

**Curently working on:**
- better projection handling
- refactoring raster class hierarchy

**Thinking about:**
- polygon holes
- memory boundedness

##CONTENTS

- vector
    - gpx_parser :  Parser for GPX files exported from GPS devices
    - geojson :     Classes and functions for handling GeoJSON files
    - vtk :         XML-based VTK interface
    - guppy :       Vector geometry classes
    - shapefile :   Snapshot of pyshp module for reading and writing shapefiles
    - shp_funcs :   Shapefile-to-guppy conversions
    - stats :       Geostatistical functions
    - xyfile :      ASCII table functions

- raster
    - grid :        Basic Grid types, including StructuredGrid and RegularGrid
    - aaigrid :     Class for reading, writing, and manipulating ESRI ASCII grids
    - flow :        Stream flow functions
    - raster :      General purpose raster math
    - streamline :  Streamline calculation
    - cfuncs :      Non-user-facing Cython module for performance

- tests : unit tests


##FORMATS

*Karta* attempts to provide a basic working interface to several of common file
formats. Currently partially supported are:

- vector
    - ASCII tables (XYZ) (r,w)
    - GeoJSON (r,w)
    - VTK (w)
    - ESRI Shapefiles (r,w)
- raster
    - ESRI ASCII Grids (r,w)

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

