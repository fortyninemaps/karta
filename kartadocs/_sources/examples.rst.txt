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


