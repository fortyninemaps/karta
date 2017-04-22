Karta Documentation
===================

*Karta* is a package for spatial analysis in Python. It streamlines data
processing by providing generic geographical types for vector and raster sources
as well as a selection of analysis functions.

*Karta* is divided into two major sections, ``karta.raster`` for gridded
image-like data, and ``karta.vector`` for point-wise geometry data. Users of GIS
software will find the distinction familiar. Rasters are typically used for
applications such as satellite imagery, elevation and hillshade maps, and
heatmaps. Vector datasets are used when it makes sense to talk about discrete
geometries, as in a network of stations, a road system, or maps administrative
boundaries.

In addition to the above, *Karta* contains the modules ``karta.crs`` and
``karta.geodesy``, which together provide support classes and functions that
help represent raster and vector object as positions on earth. There are
numerous ways of representing positions, including geographical
longitude/latitude, n-vectors, and easting/northings on various projected
systems. The CRS class encodes the type of coordinate reference system used by
an individual or set of objects, and provides a translation layer permitting
intercomparison of objects using different systems.

*Karta* is tested with Python 2.7 and Python 3.4+. Suggestions, bug reports, and
pull requests are welcome via the `Github page`_.

.. _`Github page`: https://www.github.com/fortyninemaps/karta

Getting Started
---------------

For a whirlwind tour of some basic usage of *Karta*, look through the `tutorial`_.

.. _`tutorial`: _static/tutorial.html

Data Formats
------------

*Karta* provides a basic working interface to several common file formats.
Currently supported include:

-  vector

   -  GeoJSON (r,w)
   -  ESRI shapefiles (GDAL) (r,w)
   -  GPS eXchange (GPX) (r,w)

-  raster

   -  ESRI ASCII grid (r,w)
   -  GeoTiff (GDAL) (r,w)

Furthermore, *Karta* vector geometries implement the Python geo_interface_
protocol, permitting direct data interchange with other packages that implement
it, such as Shapely_ and ArcPy_.

.. _geo_interface: https://gist.github.com/sgillies/2217756
.. _Shapely: https://github.com/Toblerity/Shapely
.. _ArcPy: http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//000v00000153000000


Examples
--------

Read or create vector geometries:

::

    point = Point((-130.0, 52.0), crs=LonLatWGS84)
    line = read_geojson("linedata.json")
    polygon = Polygon([(-515005.78, -1301130.53),
                       (-579174.89, -1282271.94),
                       (-542977.83, -1221147.82),
                       (-437864.05, -1251641.55),
                       (-438160.72, -1252421.48),
                       (-437961.28, -1285314.00)],
                       crs=NSIDCNorth)

Perform simple queries:

::

    point2 = Point((-25.0, 48.0), crs=LonLatWGS84)
    point.distance(point2)          # Distance in geographical units
    line.intersects(polygon)        # True or False
    ch = polygon.convex_hull()      # Returns a new polygon
    ch.to_shapefile("poly.shp")

Load and manipulate raster data:

::

    grid = read_gtiff("landsat_scene.tif")  # Leverages GDAL
    grid.profile(line)              # Collect data along a line
    grid.resample(500.0, 500.0)     # Return a grid resampled at a new resolution


License
-------

This software is provided under the MIT license.

MIT License:
++++++++++++

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

