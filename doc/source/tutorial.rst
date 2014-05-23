.. _tutorial:
Karta Tutorial
==============

Introduction
------------

Karta provides a lightweight set of tools for performing analyses of
geographical data. The organization of karta is around a set of
container classes for vector and raster data with builtin methods for
common tasks. This tutorial provides a brief introduction to some of the
main parts of karta.

Should you come across any mistakes, I would very much appreciate if you
could let me know, or even better, provide a pull request on
`Github <https://github.com/njwilson23/karta>`__!

Definitions
-----------

-  **vector data** refers to data can is treated as a set of connected
   or disconnected vertices. Examples might be road networks, a set of
   borders, geophysical survey lines, or the path taken by a bottle
   floating in an ocean current. In karta, these data are classified as
   belonging to *Point*, *Multipoint*, *Line* or *Polygon* classes. Some
   questions that might be asked of vector data include
-  which of these points are contained in this polygon?
-  how many times and where do these lines intersect each other?
-  what is the average distance travelled by a particle?

-  **raster data**, in contrast, is data that are typically thought of
   in terms of pixels or a grid of values covering a surface. Examples
   might be an elevation map, satellite image, or an upstream area map.
   Operations on raster data might include slope, aspect, and curvature
   calculations, up and downsampling, and interpolation.

-  **Coordinate reference system** refers to a system of relating
   measurements on a coordinate system to actual positions in space,
   e.g. on the curved surface of the Earth. karta includes very basic
   support of projected and geographical coordinates, but extending this
   system through *pyproj* is something I would like to accomplish in
   the future.

Vector data
-----------

Let's experiment with some vector data.

::

    import karta.vector as kv

The ``Point``, ``Multipoint``, ``Line``, and ``Polygon`` classes can all
be instantiated by providing vertices, and optionally, associated data
and metadata.

::

    from karta.vector import Point, Multipoint, Line, Polygon

    pt = Point((-123.1, 49.25))

    mp = Multipoint([(-122.93, 48.62),
                     (-123.10, 48.54),
                     (-122.90, 48.49),
                     (-122.81, 48.56)])

    line = Line([(-124.35713, 49.31437),
                 (-124.37857, 49.31720),
                 (-124.39442, 49.31833),
                 (-124.40311, 49.31942),
                 (-124.41052, 49.32203),
                 (-124.41681, 49.32477),
                 (-124.42278, 49.32588)])

    poly = Polygon([(-25.41, 67.03),
                    (-24.83, 62.92),
                    (-12.76, 63.15),
                    (-11.44, 66.82)])

Each geometrical object now contains a vertex/vertices in a cartesian
plane.

We may be interested in determining whether our point is within our
polygon

::

    poly.contains(pt)

    >> False

or whether our line crosses the polygon

::

    line.intersects(poly)

    >> False

The vertices of multiple vertex objects can be iterated through and
sliced:

::

    subline = line[2:-2]
    for pt_ in subline:
        print(pt_.vertex)

    >> (-124.39442, 49.31833)
       (-124.40311, 49.31942)
       (-124.41052, 49.32203)

We could have instantiated the objects by passing the ``crs`` keyword
argument a CRS instance, as in

::

    pt = Point((-123.1, 49.25), crs=kv.LONLAT)

to indicate geographical coordinates. Having done do, we could calculate
the great circle distance from our point to another point:

::

    pt2 = Point((-70.66, 41.52), crs=kv.LONLAT)
    pt.greatcircle(pt2)

    >> 4109559.587727985

**TODO: describe other geometries**

Associated data
~~~~~~~~~~~~~~~

By using the ``data`` keyword argument, additional data can be
associated with a vector geometry.

::

    mp = Multipoint([(1, 1), (3, 1), (4, 3), (2, 2)],
                    data={"species": ["T. officianale", "C. tectorum",
                                      "M. alba", "V. cracca"]})

The data can be a list or a dictionary of lists, and are propogated
through subsequent operations.

::

    pt = mp[2]
    pt

    >> Point((4, 3))

    pt.data["species"]

    >> 'M. alba'

Metadata at the geometry level rather than the point level can be
provided using the ``properties`` keyword argument, which accepts a
dictionary.

::

    poly = Polygon([(-25.41, 67.03),
                    (-24.83, 62.92),
                    (-12.76, 63.15),
                    (-11.44, 66.82)],
                   properties={"geology": "volcanic",
                               "alcohol": "brennivin"})

Visualizing and importing/exporting data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``get_coordinate_lists`` method provides lists of coordinates, which
make is easy to visualize a geometry.

::

    import matplotlib.pyplot as plt
    plt.plot(*line.get_coordinate_lists())

Data can be read from several common formats, including ESRI shapefiles
(through bindings to the *pyshp* module), GeoJSON, GPX, comma separated
value tables. Convenience functions are kept in the
``karta.vector.read`` namespace.

Saving data happens as

::

    mp.to_vtk("my_vtk.vtk")

    line.to_shapefile("my_shapefile")

    poly.to_geojson("my_json.json")

*Note: Support for various file formats is of varying quality, and I
tend to rely on GeoJSON the most.*

Raster
------

