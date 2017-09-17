Table of Contents
=================

.. do not remove the toctree directive

.. toctree::
   :maxdepth: 2

   installation
   vector
   raster
   crs
   examples

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



Indices and tables
------------------

* :ref:`genindex`

.. * :ref:`search`

