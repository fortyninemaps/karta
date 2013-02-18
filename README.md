README
------

*geo_tools* contains a collection of loosely-related Python modules for
performing lightweight geospatial data analysis. *geo_tools* replicates the
functionality of GDAL, OGR, and GSTAT in some cases, but does so in a minimal
package with few dependencies, and attempts to solve problems in ways that
would be considered 'pythonic'. Although these other packages may have
advantages in performance or capability, *geo_tools* is designed to be portable
require as little extra software infrastructure as possible.

*geo_tools* should be considered a work in progress. No assurances are provided
about results derived using it.

##CONTENTS

- vector
    - gpx_parser      Parser for GPX files exported from GPS devices
    - geojson         Classes and functions for handling GeoJSON files
    - guppy           Vector geometry classes
    - shapefile       Snapshot of pyshp module for reading and writing
    - shapefiles
    - stats           Geostatistical functions
    - xy              XY(Z) table functions

- raster
    - aaigrid_driver  Class for reading and writing ESRI ASCII Grids
    - flow            Stream flow functions
    - raster          General purpose raster math
    - streamline      Streamline calculation
    - cfuncs          Cython-module for performance, not intended to be user-facing

- tests
    - miscellaneous test scripts


##FORMATS

*geo_tools* attempts to provide a basic working interface to several of common
file formats. Currently supported are:

- vector
    - ASCII tables
    - GeoJSON
    - VTK
- raster
    - ESRI ASCII Grids


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

