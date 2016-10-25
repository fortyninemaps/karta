Package reference
=================

This is a nearly-exhaustive listing of *karta* classes and function. Some
experimental, deprecated, or private components are excluded.

Raster package (``karta.raster``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``karta.raster`` package provides *Grid* classes as well as readers for ESRI
ASCII grids and GeoTiff files. The most commonly used *Grid* subclass is
*RegularGrid*, which can efficiently handle gridded datasets with constant
offsets between nodes. *Grids* support geographically-aware clipping, sampling,
and profile extraction.

A more specialized grid class, *WarpedGrid*, is experimental.

.. automodule:: karta.raster.grid

.. autoclass:: karta.raster.grid.Grid
    :members:

RegularGrid
-----------

.. autoclass:: karta.raster.grid.RegularGrid
    :members:

WarpedGrid
----------

.. autoclass:: karta.raster.WarpedGrid
    :members:

SimpleBand
----------

.. autoclass:: karta.raster.band.SimpleBand
    :members:

CompressedBand
--------------

.. autoclass:: karta.raster.band.CompressedBand
    :members:

Miscellaneous raster functions
------------------------------

.. automodule:: karta.raster.misc
    :members:

Raster IO modules
-----------------

.. automodule:: karta.raster.read
    :members:

GeoTiff (GDAL interface)
++++++++++++++++++++++++

.. automodule:: karta.raster._gdal
    :members:

ESRI ASCII grids
++++++++++++++++

.. automodule:: karta.raster._aai
    :members:

Vector package (``karta.vector``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The package ``karta.vector`` provides a *Geometry* class, subclasses for
*Point*, *Line*, and *Polygon* types and their multipart counterparts, as well
as readers and writers for ESRI shapefiles, GeoJSON, and GPX files.

All concrete geometries support the `__geo_interface__`_ attribute, and map to
*Point*, *MultiPoint*, *LineString*, and *Polygon* types, respectively.

.. _`__geo_interface__`: https://gist.github.com/sgillies/2217756

Geometry
--------

.. automodule:: karta.vector.geometry

Geometry
++++++++

.. autoclass:: karta.vector.geometry.Geometry
    :members:

Point
+++++

.. autoclass:: karta.vector.Point
    :members:
    :inherited-members:

Line
++++

.. autoclass:: karta.vector.Line
    :members:
    :inherited-members:

Polygon
+++++++

.. autoclass:: karta.vector.Polygon
    :members:
    :inherited-members:

Multipoint
++++++++++

.. autoclass:: karta.vector.Multipoint
    :members:
    :inherited-members:

Multiline
+++++++++

.. autoclass:: karta.vector.Multiline
    :members:
    :inherited-members:

Multipolygon
++++++++++++

.. autoclass:: karta.vector.Multipolygon
    :members:
    :inherited-members:

Metadata tables
---------------

.. automodule:: karta.vector.table

.. autoclass:: karta.vector.table.Table
    :members:

.. autoclass:: karta.vector.table.Indexer
    :members:

Trees
-----

Quadtree
++++++++

.. automodule:: karta.vector.quadtree
    :members:

R-Tree
++++++

.. automodule:: karta.vector.rtree
    :members:

Vector IO modules
-----------------

.. automodule:: karta.vector.read
    :members:

GeoJSON
+++++++

.. automodule:: karta.vector.geojson
    :members:

ESRI shapefile (GDAL interface)
+++++++++++++++++++++++++++++++

.. automodule:: karta.vector.shp
    :members:

GPS Exchange
++++++++++++

.. automodule:: karta.vector.gpx
    :members:

Managing coordinate systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Coordinate reference systems describe how the numerical coordinates of a
geometry or grid relate to space in the real world. Difference coordinate
systems take different approaches to reducing position on the curved surface of
the earth to two-dimensional coordinates, making different assumptions about
Earth's geometry and choosing to preserve different invariants.

In *Karta*, coordinate systems objects handle both projection from geographical
coordinates and geodetic operations such as calculating distances and azimuths.
Basic coordinate systems (*Cartesian* and *SphericalEarth*) are provided, while
more specialized systems employ `pyproj`_ to handle calculations.

.. _pyproj: http://jswhit.github.io/pyproj/

The following predefined coordinate systems instances are available:

    - *SphericalEarth*
    - *LonLatWGS84*
    - *LonLatNAD27*
    - *LonLatNAD83*
    - *UPSNorth*
    - *UPSSouth*
    - *NSIDCNorth*
    - *NSIDCSouth*
    - *LambertEqualArea*
    - *GallPetersEqualArea*
    - *WebMercator*

Additionally, the *Cartesian* class can be used without instantiation.

To use a pre-defined coordinate system,

::

    from karta.crs import LonLatWGS84

The object ``LonLatWGS84`` can be passed to *Geometry* and *Grid* instances in
*Karta*. To define their coordinate system.

To create a specialized coordinate system, provide a *proj.4* string and a
spheroid definition.

::

    from karta.crs import Proj4CRS
    my_crs = Proj4CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 "
                      "+lon_0=-126 +x_0=1000000 +y_0=0 "
                      "+ellps=GRS80 +datum=NAD83 +units=m +no_defs",
                      "+ellps=GRS80")

See the proj.4_ documentation for more details.

.. _proj.4: http://trac.osgeo.org/proj/

CRS (``karta.crs``)
-------------------

.. automodule:: karta.crs
    :members:

Geodesy (``karta.geodesy``)
---------------------------

Functions from ``karta.geodesy`` should generally be accessed through a ``CRS``
subclass rather than directly.

.. automodule:: karta.geodesy
    :members:
