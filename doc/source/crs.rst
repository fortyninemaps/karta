Coordinate reference systems
============================

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

    from karta.crs import ProjectedCRS
    my_crs = ProjectedCRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 "
                          "+lon_0=-126 +x_0=1000000 +y_0=0 "
                          "+ellps=GRS80 +datum=NAD83 +units=m +no_defs")

See the proj.4_ documentation for more details.

.. _proj.4: http://trac.osgeo.org/proj/


