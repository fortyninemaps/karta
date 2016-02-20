# Highlights

## changes with 0.6

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
