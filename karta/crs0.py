""" Creates a database of CRS IDs. The CRSRegister class tracks a dictionary of
CRS instances. When a new CRS is requested, it is returned from the dictionary
or created on the fly.

Customized reference systems may be added at runtime as follows:

::

    from karta.crs import crsreg

    crsreg.add_CRS("my_CRS", proj=[projection parameters],
                             geod=[geod parameters],
                             crstype=[one of 'projected', 'local', or 'geographical'],
                             id=[OGC URN (string)])

For example, to add a CRS for handling data in the European Terrestrial
Reference System of 1989, use:

::

    crsreg.add_CRS("ETRS89", proj={"proj": "lonlat", "towgs84": [0,0,0,0,0,0,0]},
                             geod={"ellps": "GRS80"},
                             crstype="geographical",
                             id="urn:ogc:def:crs:EPSG::4258")
"""

import pyproj

class CRS(object):

    """ The `CRS` class represents a specific coordinate reference systems.
    Each `CRS` instance has a `crstype`, which is one of:

        - "geographical"
        - "projected"
        - "local"

    and an `id` string, which is intended to be an Open Geospatial Consortium
    uniform resource name ('uri'). Also, a `CRS` may have a callable `proj`
    method that projects between the `CRS` and geographical longitude and
    latitude, and a callable `geod` attribute that computes geodetic distances.
    """
    def __init__(self, proj=None, geod=None, crstype="unkown", id=None):
        """'proj' : pyproj.Proj or None

        'geod' : pyproj.Geod or None

        'crstype' : string
        Must be one of 'geographical', 'projected', or 'local'

        'id' : string
        Open Geospatial Consortium uniform resource name
        """
        if crstype in ("geographical", "projected", "local", "unknown"):
            self.crstype = crstype
        else:
            raise TypeError("CRS argument 'crstype' must be 'geographical', "
                            "'projected', 'local', or 'unknown'")
        self.proj = proj
        self.geod = geod
        self.id = id
        return

    def __repr__(self):
        return "{0} CRS (urn: {1})".format(self.crstype, self.id)

    def __str__(self):
        return "{0} CRS (ID: {1}".format(self.crstype, self.id)

    def __eq__(self, other):
        if (self.crstype is not "unknown" and self.crstype == other.crstype and
            self.id == other.id and
            getattr(self.proj, "srs", 0) == getattr(other.proj, "srs", 0) and
            getattr(self.geod, "initstring", 0) == getattr(other.geod, "initstring", 0)):
                return True
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

class CRSRegister(object):
    """ The `CRSRegister` is a database of known and initialized coordinate
    systems. Individual `CRS` instances can be accessed as attributes of the
    `CRSRegister`. If they are already in use, a reference to the existing
    `CRS` is returned to avoid duplication of thousands of `CRS` instances.
    Otherwise, known references systems are initialized and returned.
    """
    
    def __init__(self):
        self.register = {}

    def __getattr__(self, name):
        return self.get_CRS(name)

    def add_CRS(self, name, **kwargs):
        """ Add a ``crs.CRS`` instance with `**kwargs` with the identifier
        `name`. """
        self.register[name] = CRS(**kwargs)

    def get_CRS(self, name):
        """ Get a reference to a ``crs.CRS`` with identifier `name`.
        """
        crs = self.register.get(name, None)
        if crs is None:
            if name in CRSDict:
                cfg = CRSDict[name]

                if cfg["proj"] is None:
                    p = None
                else:
                    p = pyproj.Proj(**cfg["proj"])

                if cfg["geod"] is None:
                    g = None
                else:
                    g = pyproj.Geod(**cfg["geod"])

                crs = CRS(proj=p, geod=g, crstype=cfg["crstype"], id=cfg["id"])
                self.register[name] = crs
            else:
                raise CRSLookupError("CRS '{0}' not defined".format(name))
        return crs

CRSDict = {"CARTESIAN" : {"proj": None,
                          "geod": None,
                          "crstype": "local",
                          "id": {"urn": "urn:ogc:def:crs:EPSG::5806"}},

           # Geographical
           "LONLAT_WGS84" : {"proj": {"proj": "lonlat"},
                             "geod": {"ellps":"WGS84"},
                             "crstype": "geographical",
                             "id": {"urn": "urn:ogc:def:crs:EPSG::4326"}},
           "LONLAT_NAD27" : {"proj": {"proj": "lonlat"},
                             "geod": {"ellps":"clrk66"},
                             "crstype": "geographical",
                             "id": {"urn": "urn:ogc:def:crs:EPSG::4267"}},
           "LONLAT_NAD83" : {"proj": {"proj": "lonlat"},
                             "geod": {"ellps":"GRS80"},
                             "crstype": "geographical",
                             "id": {"urn": "urn:ogc:def:crs:EPSG::4269"}},

           # Polar stereographic
           "UPSNORTH_WGS84" : {"proj": {"proj": "stere",
                                        "lat_0": 90,
                                        "lat_ts": 90,
                                        "lon_0": 0,
                                        "k": 0.994,
                                        "x_0": 2000000,
                                        "y_0": 2000000,
                                        "units": "m",
                                        "no_defs": True},
                               "geod": {"ellps": "WGS84"},
                               "crstype": "projected",
                               "id": {"urn": "urn:ogc:def:crs:EPSG::32661"}},
           "UPSSOUTH_WGS84" : {"proj": {"proj": "stere",
                                        "lat_0": -90,
                                        "lat_ts": -90,
                                        "lon_0": 0,
                                        "k": 0.994,
                                        "x_0": 2000000,
                                        "y_0": 2000000,
                                        "units": "m",
                                        "no_defs": True},
                               "geod": {"ellps": "WGS84"},
                               "crstype": "projected",
                   "id": {"urn": "urn:ogc:def:crs:EPSG::32761"}},
           "NSIDCNORTH_WGS84" : {"proj": {"proj": "stere",
                                          "lat_0": 90,
                                          "lat_ts": 70,
                                          "lon_0": -45,
                                          "k": 1,
                                          "x_0": 0,
                                          "y_0": 0,
                                          "units": "m",
                                          "no_defs": True},
                                 "geod": {"ellps": "WGS84"},
                                 "crstype": "projected",
                                 "id": {"urn": "urn:ogc:def:crs:EPSG::3413"}},
           "NSIDCSOUTH_WGS84" : {"proj": {"proj": "stere",
                                          "lat_0": -90,
                                          "lat_ts": -70,
                                          "lon_0": 0,
                                          "k": 1,
                                          "x_0": 0,
                                          "y_0": 0,
                                          "units": "m",
                                          "no_defs": True},
                                 "geod": {"ellps": "WGS84"},
                                 "crstype": "projected",
                                 "id": {"urn": "urn:ogc:def:crs:EPSG::3976"}},

           # Unknown
           "UNKNOWN" : {"proj": None,
                        "geod": None,
                        "crstype": "local",
                        "id": {}}
           }

# Aliases
CRSDict["LONLAT"] = CRSDict["LONLAT_WGS84"]
CRSDict["UPSNORTH"] = CRSDict["UPSNORTH_WGS84"]
CRSDict["UPSSOUTH"] = CRSDict["UPSSOUTH_WGS84"]
CRSDict["NSIDCNORTH"] = CRSDict["NSIDCNORTH_WGS84"]
CRSDict["NSIDCSOUTH"] = CRSDict["NSIDCSOUTH_WGS84"]

# This is the CRS register used globally within karta
crsreg = CRSRegister()

# UTM zones
#UTM1N = CRS("UTM Zone 1 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32601")
#UTM1S = CRS("UTM Zone 1 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32701")
#UTM2N = CRS("UTM Zone 2 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32602")
#UTM2S = CRS("UTM Zone 2 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32702")
#UTM3N = CRS("UTM Zone 3 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32603")
#UTM3S = CRS("UTM Zone 3 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32703")
#UTM4N = CRS("UTM Zone 4 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32604")
#UTM4S = CRS("UTM Zone 4 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32704")
#UTM5N = CRS("UTM Zone 5 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32605")
#UTM5S = CRS("UTM Zone 5 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32705")
#UTM6N = CRS("UTM Zone 6 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32606")
#UTM6S = CRS("UTM Zone 6 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32706")
#UTM7N = CRS("UTM Zone 7 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32607")
#UTM7S = CRS("UTM Zone 7 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32707")
#UTM8N = CRS("UTM Zone 8 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32608")
#UTM8S = CRS("UTM Zone 8 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32708")
#UTM9N = CRS("UTM Zone 9 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32609")
#UTM9S = CRS("UTM Zone 9 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32709")
#UTM10N = CRS("UTM Zone 10 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32610")
#UTM10S = CRS("UTM Zone 10 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32710")
#UTM11N = CRS("UTM Zone 11 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32611")
#UTM11S = CRS("UTM Zone 11 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32711")
#UTM12N = CRS("UTM Zone 12 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32612")
#UTM12S = CRS("UTM Zone 12 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32712")
#UTM13N = CRS("UTM Zone 13 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32613")
#UTM13S = CRS("UTM Zone 13 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32713")
#UTM14N = CRS("UTM Zone 14 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32614")
#UTM14S = CRS("UTM Zone 14 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32714")
#UTM15N = CRS("UTM Zone 15 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32615")
#UTM15S = CRS("UTM Zone 15 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32715")
#UTM16N = CRS("UTM Zone 16 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32616")
#UTM16S = CRS("UTM Zone 16 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32716")
#UTM17N = CRS("UTM Zone 17 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32617")
#UTM17S = CRS("UTM Zone 17 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32717")
#UTM18N = CRS("UTM Zone 18 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32618")
#UTM18S = CRS("UTM Zone 18 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32718")
#UTM19N = CRS("UTM Zone 19 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32619")
#UTM19S = CRS("UTM Zone 19 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32719")
#UTM20N = CRS("UTM Zone 20 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32620")
#UTM20S = CRS("UTM Zone 20 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32720")
#UTM21N = CRS("UTM Zone 21 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32621")
#UTM21S = CRS("UTM Zone 21 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32721")
#UTM22N = CRS("UTM Zone 22 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32622")
#UTM22S = CRS("UTM Zone 22 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32722")
#UTM23N = CRS("UTM Zone 23 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32623")
#UTM23S = CRS("UTM Zone 23 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32723")
#UTM24N = CRS("UTM Zone 24 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32624")
#UTM24S = CRS("UTM Zone 24 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32724")
#UTM25N = CRS("UTM Zone 25 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32625")
#UTM25S = CRS("UTM Zone 25 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32725")
#UTM26N = CRS("UTM Zone 26 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32626")
#UTM26S = CRS("UTM Zone 26 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32726")
#UTM27N = CRS("UTM Zone 27 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32627")
#UTM27S = CRS("UTM Zone 27 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32727")
#UTM28N = CRS("UTM Zone 28 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32628")
#UTM28S = CRS("UTM Zone 28 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32728")
#UTM29N = CRS("UTM Zone 29 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32629")
#UTM29S = CRS("UTM Zone 29 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32729")
#UTM30N = CRS("UTM Zone 30 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32630")
#UTM30S = CRS("UTM Zone 30 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32730")
#UTM31N = CRS("UTM Zone 31 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32631")
#UTM31S = CRS("UTM Zone 31 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32731")
#UTM32N = CRS("UTM Zone 32 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32632")
#UTM32S = CRS("UTM Zone 32 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32732")
#UTM33N = CRS("UTM Zone 33 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32633")
#UTM33S = CRS("UTM Zone 33 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32733")
#UTM34N = CRS("UTM Zone 34 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32634")
#UTM34S = CRS("UTM Zone 34 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32734")
#UTM35N = CRS("UTM Zone 35 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32635")
#UTM35S = CRS("UTM Zone 35 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32735")
#UTM36N = CRS("UTM Zone 36 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32636")
#UTM36S = CRS("UTM Zone 36 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32736")
#UTM37N = CRS("UTM Zone 37 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32637")
#UTM37S = CRS("UTM Zone 37 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32737")
#UTM38N = CRS("UTM Zone 38 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32638")
#UTM38S = CRS("UTM Zone 38 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32738")
#UTM39N = CRS("UTM Zone 39 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32639")
#UTM39S = CRS("UTM Zone 39 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32739")
#UTM40N = CRS("UTM Zone 40 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32640")
#UTM40S = CRS("UTM Zone 40 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32740")
#UTM41N = CRS("UTM Zone 41 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32641")
#UTM41S = CRS("UTM Zone 41 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32741")
#UTM42N = CRS("UTM Zone 42 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32642")
#UTM42S = CRS("UTM Zone 42 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32742")
#UTM43N = CRS("UTM Zone 43 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32643")
#UTM43S = CRS("UTM Zone 43 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32743")
#UTM44N = CRS("UTM Zone 44 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32644")
#UTM44S = CRS("UTM Zone 44 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32744")
#UTM45N = CRS("UTM Zone 45 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32645")
#UTM45S = CRS("UTM Zone 45 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32745")
#UTM46N = CRS("UTM Zone 46 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32646")
#UTM46S = CRS("UTM Zone 46 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32746")
#UTM47N = CRS("UTM Zone 47 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32647")
#UTM47S = CRS("UTM Zone 47 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32747")
#UTM48N = CRS("UTM Zone 48 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32648")
#UTM48S = CRS("UTM Zone 48 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32748")
#UTM49N = CRS("UTM Zone 49 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32649")
#UTM49S = CRS("UTM Zone 49 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32749")
#UTM50N = CRS("UTM Zone 50 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32650")
#UTM50S = CRS("UTM Zone 50 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32750")
#UTM51N = CRS("UTM Zone 51 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32651")
#UTM51S = CRS("UTM Zone 51 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32751")
#UTM52N = CRS("UTM Zone 52 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32652")
#UTM52S = CRS("UTM Zone 52 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32752")
#UTM53N = CRS("UTM Zone 53 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32653")
#UTM53S = CRS("UTM Zone 53 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32753")
#UTM54N = CRS("UTM Zone 54 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32654")
#UTM54S = CRS("UTM Zone 54 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32754")
#UTM55N = CRS("UTM Zone 55 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32655")
#UTM55S = CRS("UTM Zone 55 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32755")
#UTM56N = CRS("UTM Zone 56 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32656")
#UTM56S = CRS("UTM Zone 56 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32756")
#UTM57N = CRS("UTM Zone 57 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32657")
#UTM57S = CRS("UTM Zone 57 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32757")
#UTM58N = CRS("UTM Zone 58 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32658")
#UTM58S = CRS("UTM Zone 58 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32758")
#UTM59N = CRS("UTM Zone 59 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32659")
#UTM59S = CRS("UTM Zone 59 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32759")
#UTM60N = CRS("UTM Zone 60 North (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32660")
#UTM60S = CRS("UTM Zone 60 South (WGS 84)", "projected", "urn:ogc:def:crs:EPSG::32760")

class CRSError(Exception):
    """ Exception to raise for invalid geodetic operations. """
    def __init__(self, message=''):
        self.message = message

 
class CRSLookupError(Exception):
    def __init__(self, message=''):
        self.message = message

