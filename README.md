[![Build Status](https://travis-ci.org/fortyninemaps/karta.svg?branch=master)](https://travis-ci.org/fortyninemaps/karta)
[![Build status](https://ci.appveyor.com/api/projects/status/viiimwp5pu7ff2bp?svg=true)](https://ci.appveyor.com/project/njwilson23/karta)
[![Coverage Status](https://coveralls.io/repos/github/fortyninemaps/karta/badge.svg?branch=master)](https://coveralls.io/github/fortyninemaps/karta?branch=master)

![Karta](https://raw.githubusercontent.com/fortyninemaps/karta/gh-pages/images/karta_logo.png)

*Karta* is a package for spatial analysis in Python. It simplifies geospatial
data processing by providing efficient generic classes for vector and raster
data sources, as well as a selection of analysis functions.

## Vector data types

Data are represented as `Point`, `Line`, `Polygon`, `Multipoint`, `Multiline`,
and `Multipolygon` instances.

All data contain a `.crs` member encoding coordinate reference information. All
vector geometries possess a `.properties` dict containing free-form metadata.
*Multipart* geometries additionally possess a `.data` member which is a simple
typed table-like data structure.

Geometries implement methods for computing distances, directions, and spatial
characteristics. *Multipart* geometries support fast spatial indexing and
queries.

GeoJSON and ESRI shapefile formats are supported for reading and writing.
Experimental support for GPX XML files is in the `karta.vector.gpx` submodule.

Vector geometries implement the Python [`__geo_interface__`
attribute](https://gist.github.com/sgillies/2217756) for vector geometries. This
permits data to be exchanged between *Karta* and external modules that also
implement `__geo_interface__` (e.g.
[shapely](https://github.com/Toblerity/Shapely),
[fastkml](https://fastkml.readthedocs.org/en/latest/)).

## Raster data types

The primary raster data type is the `RegularGrid`, which represents one or more
2-d arrays of pixels spaced via an affine transformation. `RegularGrids` are
backed by one of several `Band` implementations, with the default implementation
using the [blosc](http://www.blosc.org/) compression library for efficient
in-memory storage. There is experimental support for disk-backed storage via
GDAL.

Grids may be queried, resampled, sliced, masked, and merged. Arbitrary
array-based functions may be mapped to raster data with `RegularGrid.apply()`.
Raster functions including slope, gradient, and hillshade are in the
`karta.raster.misc` submodule.

GeoTIFF images are the primary supported format, however ESRI ASCII grids may
also be used (with limitations due to the format).

## Coordinate reference systems

Data in Karta is referenced to positions on earth via `CRS` objects that
implement projection and geodetic methods. Coordinate reference systems may be
either geographical or projected.

**Geographical CRS** objects return spatial relationships in terms of the true
computed distances and azimuths on a spherical or ellipsoidal Earth.

**Projected CRS** objects (e.g. UTM, Polar Stereographic, and Web Mercator)
return spatial relationships in terms of a flat plane, dependent on the
projection.

## Examples

Read or create vector geometries:

```python
point = Point((-130.0, 52.0), crs=LonLatWGS84)
line = read_geojson("linedata.json")
polygon = Polygon([(-515005.78, -1301130.53),
                   (-579174.89, -1282271.94),
                   (-542977.83, -1221147.82),
                   (-437864.05, -1251641.55),
                   (-438160.72, -1252421.48),
                   (-437961.28, -1285314.00)],
                   crs=NSIDCNorth)
```
Perform simple queries:
```python
point2 = Point((-25.0, 48.0), crs=LonLatWGS84)
point.distance(point2)          # Distance in geographical units
line.intersects(polygon)        # True or False
ch = polygon.convex_hull()      # Returns a new polygon
ch.to_shapefile("poly.shp")
```
Load and manipulate raster data:
```python
grid = read_gtiff("landsat_scene.tif")  # Leverages GDAL
grid.profile(line)              # Collect data along a line
grid.resample(500.0, 500.0)     # Return a grid resampled at a new resolution
```

## Installation

*Karta* currently supports Python 2.7 and Python 3.4+.

The easiest way to install is via `pip`. Installation requires a recent version
of `setuptools`.

    pip install -U setuptools
    pip install karta

### Building from source

Building requires Cython and a C99-compliant compiler.

    pip install Cython

Clone the repository and build.

    git clone https://github.com/fortyninemaps/karta.git karta
    cd karta/
    python setup.py build

## Documentation

See the [online manual](http://karta.fortyninemaps.com/kartadocs/karta-manual.html),
the [tutorial](http://karta.fortyninemaps.com/kartadocs/_static/tutorial.html), or read the
[API documentation](http://karta.fortyninemaps.com/kartadocs/reference.html).

## Contributing

The source code lives at
[github.com/fortyninemaps/karta](github.com/fortyninemaps/karta).

Bug reports, feature requests, and pull requests are welcome.

Run unit tests with `python tests/runtests.py`.

The manual is built using [Sphinx](http://sphinx-doc.org/) and requires
[numpydoc](https://github.com/numpy/numpydoc).

