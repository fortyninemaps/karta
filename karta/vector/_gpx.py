"""
GPX IO and manipulation (UNDER DEVELOPMENT)

Overview
--------

`GPX` represents a GPX document that can be written to or read from

To do:
    - XML namespaces aren't really handled properly
    - metadata node not addressed
"""

import collections
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import ElementTree, Element
from xml.dom import minidom

Point = collections.namedtuple("Point", ["lonlat", "properties", "extensions"])
Trkseg = collections.namedtuple("Trkseg", ["trkpts", "properties", "extensions"])
Track = collections.namedtuple("Track", ["trksegs", "properties", "extensions"])
Route = collections.namedtuple("Route", ["rtepts", "properties", "extensions"])

ET.register_namespace("", "http://www.topografix.com/GPX/1/1")
ns = "{http://www.topografix.com/GPX/1/1}"

VALID_PROPERTIES = ("ele", "time", "magvar", "geoidheight", "name", "cmt",
                    "desc", "src", "link", "sym", "type", "fix", "sat", "hdop",
                    "vdop", "pdop", "ageofdgpsdata", "dgpsid")

def strip_namespace(s):
    return s[s.index("}")+1:]

class GPX(object):
    """ Represents a GPX document, with an internal representation of
    waypoints, tracks, and route that loosely approximates the XML structure.
    Provides methods to easily add Point-like and Line-like objects as GPX
    types. """

    def __init__(self, f=None, waypoints=None, tracks=None, routes=None):
        """ Create a GPX object, either from a GPX file or from lists of
        waypoints, tracks, and routes. """

        self.waypts = []
        self.tracks = []
        self.routes = []

        if f is not None:
            self.fromfile(f)

        else:
            if waypoints is not None:
                for waypt in waypoints:
                    self.add_waypoint(waypt)
            if tracks is not None:
                for track in tracks:
                    self.add_track(track)
            if routes is not None:
                for route in routes:
                    self.add_route(route)
        return

    @staticmethod
    def _readextensions(node):
        extensions = {}
        try:
            for ext in node.find(ns + "extensions"):
                extensions[strip_namespace(ext.tag)] = ext.text
        except TypeError:
            pass
        return extensions

    @staticmethod
    def _readproperties(node, exclude=()):
        properties = {}
        for subnode in node:
            tag = strip_namespace(subnode.tag)
            if tag not in exclude:
                properties[tag] = subnode.text
        return properties

    def _readwpt(self, wpt):
        properties = self._readproperties(wpt, exclude=("extensions",))
        extensions = self._readextensions(wpt)
        lon = round(float(wpt.attrib["lon"]), 6)
        lat = round(float(wpt.attrib["lat"]), 6)
        return Point((lon, lat), properties, extensions)

    @staticmethod
    def _dict2gpx(parent, properties):
        for p in properties:
            sub = Element(ns + p)
            sub.text = str(properties[p])
            parent.append(sub)
        return parent

    def _extensions2gpx(self, parent, extensions):
        ext = Element(ns + "extensions")
        ext = self._dict2gpx(ext, extensions)
        parent.append(ext)
        return parent

    def _build_gpx_wpt(self, waypt, tag="wpt"):
        """ Build <wpt> node. """
        wpt = Element(ns + tag, lon=str(waypt.lonlat[0]), lat=str(waypt.lonlat[1]))
        wpt = self._dict2gpx(wpt, waypt.properties)
        wpt = self._extensions2gpx(wpt, waypt.extensions)
        return wpt

    def _build_gpx_trk(self, track):
        """ Build "trk" nodes. """
        trk = Element(ns + "trk")
        trk = self._dict2gpx(trk, track.properties)
        trk = self._extensions2gpx(trk, track.extensions)

        for segment in track.trksegs:
            trkseg = Element(ns + "trkseg")
            trkseg = self._dict2gpx(trkseg, segment.properties)
            trkseg = self._extensions2gpx(trkseg, segment.extensions)

            for trackpt in segment.trkpts:
                trkpt = self._build_gpx_wpt(trackpt, tag="trkpt")
                trkseg.append(trkpt)

            trk.append(trkseg)
        return trk

    def _build_gpx_rte(self, route):
        rte = Element(ns + "rte")
        rte = self._dict2gpx(rte, route.properties)
        rte = self._extensions2gpx(rte, route.extensions)
        for routept in route.rtepts:
            rte.append(self._build_gpx_wpt(routept, tag="rtept"))
        return rte

    def fromfile(self, f):
        """ Read a GPX document from *f*, which may be a filename or a
        file-like object. """

        gpxtree = ElementTree(file=f)
        self.gpx = gpxtree.getroot()

        for node in self.gpx.findall(ns + "wpt"):
            self.parse_wpt(node)

        for node in self.gpx.findall(ns + "trk"):
            self.parse_trk(node)

        for node in self.gpx.findall(ns + "rte"):
            self.parse_rte(node)
        return

    def parse_wpt(self, wpt):
        """ Parse a <wpt> node, updating self.waypoints. """
        point = self._readwpt(wpt)
        name = wpt.properties.get("name", "waypoint_" + str(len(self.waypoints)))
        self.waypoints[name] = point
        return

    def parse_trk(self, trk):
        """ Parse a <trk> node, updating self.tracks. """
        segments = []
        trkproperties = self._readproperties(trk, exclude=("trkseg",))
        trkextensions = self._readextensions(trk)

        for trkseg in trk.findall(ns + "trkseg"):

            points = [self._readwpt(trkpt) for trkpt in trkseg.findall(ns + "trkpt")]
            properties = self._readproperties(trkseg, exclude=("trkpt",))
            extensions = self._readextensions(trkseg)
            segments.append(Trkseg(points, properties, extensions))

        self.tracks.append(Track(segments, trkproperties, trkextensions))
        return

    def parse_rte(self, rte):
        properties = self._readproperties(rte, exclude=("rtept",))
        extensions = self._readextensions(rte)
        points = [self._readwpt(rtept) for rtept in rte.findall(ns + "rtept")]
        self.routes.append(Route(points, properties, extensions))
        return

    def add_waypoint(self, waypoint):
        """ Add a Point-like object as a waypoint. Properties and extension
        types are taken from waypoint.properties attribute. """
        properties = {}
        extensions = {}
        for key in waypoint.properties:
            if key in VALID_PROPERTIES:
                properties[key] = str(waypoint.properties[key])
            else:
                extensions[key] = str(waypoint.properties[key])
        waypt = Point(waypoint.vertex, properties, extensions)
        self.waypts.append(waypt)
        return

    def add_track(self, *tracks, **kw):
        """ Add Multipoint-like objects as segments of a track. Dictionaries of
        properties and extension types for the track are accepted as keyword
        arguments.

        Properties and extension types for the track segments are taken from
        the `properties` attribute of each Multipoint-like object.

        Properties and extensions types for each track point are taken from the
        `data` attribute of each Multipoint-like object.

        INCOMPLETE:

        Needs to distinguish between properties and extensions at the trkseg
        and trkpt levels

        Needs to properly add and <ele> property for Lines of rank 3
        """
        segments = []
        properties = {}
        extensions = {}

        attributes = kw.get("attributes", None)
        if attributes is not None:
            for key in attributes:
                if key in VALID_PROPERTIES:
                    properties[key] = str(attributes[key])
                else:
                    properties[key] = str(attributes[key])

        for line in tracks:
            points = []
            sg_properties = {}
            sg_extensions = {}

            for key in line.properties:
                if key in VALID_PROPERTIES:
                    sg_properties[key] = str(line.properties[key])
                else:
                    sg_extensions[key] = str(line.properties[key])

            if line.data is not None:
                pt_properties = [k for k in line.data.fields if k in VALID_PROPERTIES]
                pt_extensions = [k for k in line.data.fields if k not in VALID_PROPERTIES]
            else:
                pt_properties = []
                pt_extensions = []

            for i, xy in enumerate(line.vertices):

                ptprop = {}
                ptexte = {}
                for key in pt_properties:
                    ptprop[key] = line.d[key][i]
                    if len(xy) == 3:
                        ptprop["ele"] = xy[2]
                for key in pt_extensions:
                    ptexte[key] = line.d[key][i]

                points.append(Point((xy[0], xy[1]), ptprop, ptexte))

            segments.append(Trkseg(points, sg_properties, sg_extensions))

        self.tracks.append(Track(segments, properties, extensions))
        return

    def add_route(self, route):
        """ Add a Multipoint-like object as a route. Properties and extension
        types for the route are taken from the `properties` attribute of the
        Multipoint-like object.

        Properties and extensions types for each route point are taken from the
        `data` attribute of each Line-like object.

        INCOMPLETE:

        Needs to distinguish between properties and extensions at the trkseg
        and trkpt levels

        Needs to properly add and <ele> property for Lines of rank 3
        """
        points = []
        for i, vertex in enumerate(route.vertices):
            prop = {}
            if route.data is not None:
                for k in route.data.fields:
                    prop[k] = route.data[k][i]
            points.append(Point((vertex[0], vertex[1]), prop, {}))
        self.routes.append(Route(points, route.properties, {}))
        return

    def as_string(self, waypts=True, tracks=True, routes=True):
        """ Write GPX object to a string. Writes all waypoints, tracks, and
        routes by default, which can be changed by changing the kwargs to
        False. """
        gpx = Element(ns + "gpx", version="1.1", creator="karta")

        if waypts:
            for waypt in self.waypts:
                gpx.append(self._build_gpx_wpt(waypt))
        if tracks:
            for track in self.tracks:
                gpx.append(self._build_gpx_trk(track))
        if routes:
            for route in self.routes:
                gpx.append(self._build_gpx_rte(route))

        xmlstring = ET.tostring(gpx).decode("ascii")
        return xmlstring

    def writefile(self, fnm, **kw):
        with open(fnm, "w") as f:
            f.write(self.as_string(**kw))
        return

