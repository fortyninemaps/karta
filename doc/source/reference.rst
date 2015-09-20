Package reference
=================

Coordinate reference systems (``karta.crs``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

.. automodule:: karta.crs
    :members:

Vector package (``karta.vector``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The package ``karta.vector`` provides *Geometry* classes as well as readers and
writers for ESRI shapefiles, GeoJSON, and GPX files.

Geometries may be one of *Point*, *Multipoint*, *Line*, or *Polygon*. All
geometries support the `__geo_interface__`_ attribute, and map to *Point*,
*MultiPoint*, *LineString*, and *Polygon* types, respectively.

.. _`__geo_interface__`: https://gist.github.com/sgillies/2217756

.. automodule:: karta.vector.geometry

.. autoclass:: karta.vector.geometry.Geometry
    :members:
    :inherited-members:

.. autoclass:: karta.vector.Point
    :members:
    :inherited-members:

.. autoclass:: karta.vector.Multipoint
    :members:
    :inherited-members:

.. autoclass:: karta.vector.Line
    :members:
    :inherited-members:

.. autoclass:: karta.vector.Polygon
    :members:
    :inherited-members:

.. automodule:: karta.vector.metadata
    :members:

Vector IO modules
-----------------

.. automodule:: karta.vector.read
    :members:

.. automodule:: karta.vector.geojson
    :members:

.. automodule:: karta.vector.shp
    :members:

Raster package (``karta.raster``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``karta.raster`` package provides *Grid* classes as well as readers for ESRI
ASCII grids and GeoTiff files. The most commonly used *Grid* subclass is
*RegularGrid*, which can efficiently handle gridded datasets with constant
offsets between nodes. *Grids* support geographically-aware clipping, sampling,
and profile extraction.

A more specialized grid class, *WarpedGrid*, is experimental.

.. automodule:: karta.raster.grid
    :members:

.. automodule:: karta.raster.read
    :members:

.. automodule:: karta.raster.aaigrid
    :members:

.. automodule:: karta.raster.misc
    :members:

