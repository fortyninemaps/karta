# Highlights

## changes with 0.7.2

- fix: CompressedBand now always returns array when indexed with a slice
- fix: output dimension of grid sampling methods is now always one larger than
  the input
- `merge` now works with multiband grids
- use more efficient algorithm for intersection detection
- new `read_gtiffs` function for reading multiple files as bands
- implement `RegularGrid.coordinates()` method returning a `CoordGenerator` object
- implement `get_coordinate_lists` for Multiline and Multipolygon
- performance optimizations and fixes

## changes with 0.7.1

- better support multiband arrays for `BandIndexer`
- fix to `mask_poly` function

## changes with 0.7

- major refactor of geometry hierarchy
    - adds Multiline and Multipolygon classes
    - polygon vertices are double-nested in order to represent sub-polygons
    - only multipart geometries now contain a data attribute
- GeoJSON now attempts to project all data to geographical WGS84 coordinates as
  per https://datatracker.ietf.org/doc/draft-ietf-geojson/
    - `force_wgs84=False` can be set to avoid this behaviour
- multivertex geometries now use CoordString as a data backend, which should
  permit faster operations on the C side
- rename `Line.subsection` -> `Line.to_npoints`
- new R-tree implementation in C
- `projected` keyword argument added to `Point.distance`, `Point.azimuth`, and
  `Point.walk` to control whether geodetic or planar algorithms are used with
  the CRS is not geographical
- replaced polygon membership testing algorithm with a winding number-based test
- implemented fallback algorithm for polygon membership testing in polar regions

## changes with 0.6

- implement raster bands
- redesign raster interface and introduce [:,:] array syntax rather than values
  attribute
- replace pyshp with OGR
- deprecate some CRS classes in favour of Proj4-backed CRS
- nodata improvements for raster grids

## changes with 0.5.3

- R-tree implementation for indexing multipoint geometries
- bounding box caching for performance

## changes with 0.5

- python implementation of geodetic solutions on spheres and ellipsoids
- geodetically-correct area calculations with GeographicalCRS
- more complete geojson interface supporting nested GeometryCollections
- no longer supporting pure-Python mode (C-extensions required)
- `RegularGrid.clip` performance improvements
- added `Multipoint.convex_hull()`, returning a Polygon
- CRS objects can often report WKT and proj.4 strings (mostly depends on osgeo)
- continuing to improve CRS support in GeoJSON driver
- automatic CRS support for shapefiles
- Multipoint now uses quadtree internally for certain methods

## Changes with 0.4

- renamed `vector.guppy` -> `vector.geometry`
- `NODATA` support for raster grids
- raster grid clipping
- `gtiffread` based on GDAL
- refactor `karta.crs` with a simpler interface
- affine transforms of vector geometries

## Changes with 0.3

- `GeoJSONReader` simplifications
- implement `__geo_interface__`

## Changes with 0.2 series

- quadtrees for point data
- implemented `karta.crs` to handle coordinate reference systems
- shapefile IO improvements
