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

class GPX(object):
    """ Represents a GPX documents, with waypoints, tracks, and routes as
    attributes. """

    waypts = {}
    tracks = {}
    routes = {}

    def __init__(self, f=None, waypts=None, tracks=None, routes=None):
        """ Create a GPX object, either from a GPX file or from lists of
        waypoints, tracks, and routes. """
        if f is not None:
            self.fromfile(f)

        else:
            self.gpx = Element("gpx")

            if waypts is not None:
                self.build_wpt(waypts)
            if tracks is not None:
                self.build_trk(tracks)
            if routes is not None:
                self.build_rte(routes)

        return

    def fromfile(self, f):
        """ Read a GPX document from *f*, which may be a filename or a
        file-like object. """

        gpxtree = ElementTree(file=f)
        self.gpx = gpxtree.find("gpx")

        for el in gpxtree.findall("wpt"):
            self.parse_wpt(el)

        for el in gpxtree.findall("trk"):
            self.parse_trk(el)

        for el in gpxtree.findall("rte"):
            self.parse_rte(el)

        return

    def parse_wpt(self, node):
        pass

    def parse_trk(self, node):
        """ Parse a <trk> node, updating self.tracks. """
        name = node.getElementsByTagName('name')[0].firstChild.data

        if not name in self.tracks:
            self.tracks[name] = {}

        segments = []

        for trkseg in node.getElementsByTagName('trkseg'):

            lats = []
            lons = []
            eles = []
            coords = []
            time = []

            points = []

            for trkpt in trkseg.getElementsByTagName('trkpt'):
                lat = float(trkpt.getAttribute('lat'))
                lon = float(trkpt.getAttribute('lon'))
                ele = float(trkpt.getElementsByTagName('ele')[0].firstChild.data)
                time = trkpt.getElementsByTagName('time')[0].firstChild.data
                points.append(Point(lon, lat, ele, time))

            segments.append([a for a in points])

        self.tracks[name] = segments
        return

    def parse_rte(self, node):
        pass

    def build_wpt(self, waypts):
        pass

    def build_trk(self, tracks):
        """ Build "trk" nodes. """
        for track in tracks:

            trk = Element("trk")
            if "name" in track.properties:
                name = Element("name", text=track.data["name"])
                trk.append(name)

            for segment in track.segments:
                trkseg = Element("trkseg")
                trk.append(trkseg)

                for i, c in enumerate(segment.vertices):

                    trkpt = Element("trkpt", attrib={"lon":str(c[0]),
                                                     "lat":str(c[1])})
                    if len(c) > 2:
                        ele = Element("ele", text=str(c[2]))
                        trkpt.append(ele)
                    if "time" in segment.data:
                        time = Element("time", text=str(segment.data["time"][i]))
                        trkpt.append(time)
                    for field in segment.data:
                        node = Element(field, text=str(segments.data[field][i]))
                        trkpt.append(node)

                    trkseg.append(trkpt)
            self.gpx.append(trk)
        return

    def build_rte(self, routes):
        pass

    def writefile(self, f, waypts=True, tracks=True, routes=True):
        """ Write GPX object to a GPX file. Writes all waypoints, tracks, and
        routes by default, which can be changed by changing the kwargs to
        False. """
        ElementTree(element=gpx).write(f)
        return

