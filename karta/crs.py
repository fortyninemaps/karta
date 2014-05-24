""" This is a database of CRS IDs. """

from collections import namedtuple

CRS = namedtuple("CRS", ["name", "type", "urn"])

# List of the available coordinate systems
CARTESIAN = CRS("cartesian", "local", "urn:ogc:def:crs:EPSG::5806")
LONLAT_WGS84 = CRS("lonlat", "geographical", "urn:ogc:def:crs:EPSG::4326")

UPSNORTH_WGS84 = CRS("UPS North", "projected", "urn:ogc:def:crs:EPSG::32661")
UPSSOUTH_WGS84 = CRS("UPS South", "projected", "urn:ogc:def:crs:EPSG::32761")

NSIDCNORTH = CRS("NSIDC Sea Ice Sterographic North", "projected", "urn:ogc:def:crs:EPSG::3411")
NSIDCSOUTH = CRS("NSIDC Sea Ice Sterographic South", "projected", "urn:ogc:def:crs:EPSG::3412")

# Convenience aliases
LONLAT = LONLAT_WGS84       # Convenience



localnames = locals()
crslist = [c for c in localnames.values() if isinstance(c, CRS)]

