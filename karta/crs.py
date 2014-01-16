""" This is a database of CRS IDs. """

from collections import namedtuple

CRS = namedtuple("CRS", ["id", "type"])

LONLAT = CRS("lonlat", "geographical")
CARTESIAN = CRS("cartesian", "projected")

LONLAT = "lonlat"       # This is a temporary, nonstandard SRID

