""" This is a database of CRS IDs. """

from collections import namedtuple

CRS = namedtuple("CRS", ["shortname", "type", "urn"])

# List of the available coordinate systems
LONLAT_WGS84 = CRS("WGS84", "geographical", "urn:ogc:def:crs:EPSG::4326")
CARTESIAN = CRS("cartesian", "local", "urn:ogc:def:crs:EPSG::5806")



# Convenience aliases
LONLAT = LONLAT_WGS84       # Convenience
