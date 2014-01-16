""" This is a database of CRS IDs. """

from collections import namedtuple

CRS = namedtuple("CRS", ["id", "type"])

# List of the available coordinate systems
LONLAT = CRS("lonlat", "geographical")
CARTESIAN = CRS("cartesian", "projected")

