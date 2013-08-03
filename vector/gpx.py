"""
GPX IO and manipulation (UNDER DEVELOPMENT)

This is a rewrite of gpxparser.py, designed to fit better with guppy types.
"""

import sys
from xml.dom import minidom, Node
from xml.etree.ElementTree import ElementTree, Element
import collections

Point = collections.namedtuple("Point", ["vertex", "data", "properties"])
Trkseg = collections.namedtuple("Trkseg", ["vertices", "data", "properties"])
Track = collections.namedtuple("Track", ["segments", "properties"])
Route = collections.namedtuple("Route", ["vertices", "data", "properties"])

ns = "{http://www.topografix.com/GPX/1/1}"

def strip_namespace(s):
    return s[s.index("}")+1:]

class GPX(object):
    """ Represents a GPX documents, with waypoints, tracks, and routes as
    attributes. """

    def __init__(self, f=None, waypts=None, tracks=None, routes=None):
        """ Create a GPX object, either from a GPX file or from lists of
        waypoints, tracks, and routes. """

        self.waypts = {}
        self.tracks = {}
        self.routes = {}

        if f is not None:
            self.fromfile(f)

        else:
            self.gpx = Element("gpx", attrib={"version":"1.1",
                                              "creator":"karta"})

            if waypts is not None:
                for waypt in waypts:
                    self.add_wpt(waypt)
            if tracks is not None:
                for track in tracks:
                    self.add_trk(track)
            if routes is not None:
                for route in routes:
                    self.add_rte(route)

        return

    def fromfile(self, f):
        """ Read a GPX document from *f*, which may be a filename or a
        file-like object. """

        gpxtree = ElementTree(file=f)
        self.gpx = gpxtree.getroot()

        #for node in self.gpx:
        #    if node.tag == "wpt":
        #        self.parse_wpt(node)
        #    elif node.tag == "trk":
        #        self.parse_trk(node)
        #    elif node.tag == "rte":
        #        self.parse_rte(node)

        for node in self.gpx.findall(ns+"wpt"):
            self.parse_wpt(node)

        for node in self.gpx.findall(ns+"trk"):
            self.parse_trk(node)

        for node in self.gpx.findall(ns+"rte"):
            self.parse_rte(node)

        return

    def parse_wpt(self, wpt):
        pass

    def parse_trk(self, trk):
        """ Parse a <trk> node, updating self.tracks. """
        name = trk.get("name", str(len(self.tracks)))
        segments = []

        for trkseg in trk.findall(ns+"trkseg"):

            vertices = []
            data = {}
            properties = {}

            for node in trkseg.find(ns+"trkpt"):
                tag = strip_namespace(node.tag)
                data[tag] = []

            for trkpt in trkseg.findall(ns+"trkpt"):
                vertices.append((trkpt.attrib["lon"], trkpt.attrib["lat"]))
                for node in trkpt:
                    tag = strip_namespace(node.tag)
                    if tag in data:
                        data[tag].append(node.text)
                    else:
                        raise InconsistentFieldError("child node {0} is not "
                                        "used consistently".format(node.tag))

            try:
                for node in trkseg.find("extensions"):
                    properties[node.tag] = node.text
            except TypeError:
                pass

            segment = Trkseg(vertices, data, properties)
            segments.append(segment)

        track = Track(segments, name)
        self.tracks[name] = track
        return

    def parse_rte(self, rte):
        pass

    def add_wpt(self, waypt):
        pass

    def add_trk(self, track):
        """ Build "trk" nodes. """
        trk_name = track.properties.get("name", str(len(self.tracks)))
        self.tracks[trk_name] = track

        trk = Element(ns+"trk")
        name = Element("name")
        name.text = trk_name
        trk.append(name)

        for segment in track.segments:
            trkseg = Element(ns+"trkseg")
            trk.append(trkseg)

            for i, c in enumerate(segment.vertices):

                trkpt = Element(ns+"trkpt", attrib={"lon":str(c[0]),
                                                 "lat":str(c[1])})
                if len(c) > 2:
                    ele = Element("ele")
                    ele.text = str(c[2])
                    trkpt.append(ele)
                for field in segment.data:
                    node = Element(field)
                    node.text = str(segment.data[field][i])
                    trkpt.append(node)

                trkseg.append(trkpt)

        self.gpx.append(trk)
        return

    def add_rte(self, route):
        pass

    def writefile(self, f, waypts=True, tracks=True, routes=True):
        """ Write GPX object to a GPX file. Writes all waypoints, tracks, and
        routes by default, which can be changed by changing the kwargs to
        False. """
        ElementTree(element=self.gpx).write(f)
        return

class InconsistentFieldError(Exception):
    def __init__(self, message="No message"):
        self.message = message
    def __str__(self):
        return self.message

