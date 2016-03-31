Introduction
============

*Karta* is a simple to use package for spatial analysis in Python. It
streamlines processing of vector and raster data.

*Karta* is divided into two major sections, ``karta.raster`` and
``karta.vector``, handling gridded image-like data and point-wise geometry data
respectively. Users of GIS software will find the distinction familiar. Rasters
are typically used for applications such as satellite imagery, elevation and
hillshade maps, and heatmaps. Vector datasets are used when it makes sense to
talk about individual geometries, as in a network of stations, a road system, or
maps administrative boundaries.

In addition to the above, *karta* contains the modules ``karta.crs`` and
``karta.geodesy``, which together provide support classes and functions that
help represent raster and vector object as positions on earth. There are
numerous aways of representing positions, including geographical
longitude/latitude, n-vectors, and easting/northings on various projected
systems. The CRS class encodes the type of coordinate reference system used by
an individual or set of objects, and provides a translation layer permitting
intercomparison of objects using diferent systems.

*Karta* is tested with Python 2.7 and Python 3.4+. Suggestions, bug reports, and
pull requests are welcome, and accepted via the `Github page`_.

.. _`Github page`: https://www.github.com/fortyninemaps/karta

Data Formats
------------

*Karta* provides a basic working interface to several common file formats.
Currently supported are:

-  vector

   -  GeoJSON (r,w)
   -  ESRI shapefiles (GDAL) (r,w)
   -  GPS eXchange (GPX) (r,w)

-  raster

   -  ESRI ASCII grid (r,w)
   -  GeoTiff (GDAL) (r,w)

Furthermore, *karta* vector geometries implement the Python geo_interface_
protocol, permitting direct data interchange with other packages that implement
it, such as Shapely_ and ArcPy_.

.. _geo_interface: https://gist.github.com/sgillies/2217756
.. _Shapely: https://github.com/Toblerity/Shapely
.. _ArcPy: http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//000v00000153000000

Important classes
-----------------

This is a listing of some of the more important and commonly-used classes. A
complete listing of classes, methods, and functions is available in the `package
documentation`_. 

.. _package documentation: reference.html


Examples
--------

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

