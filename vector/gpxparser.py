"""
GPX Parser

Based on code by Rui Carmo, available at
http://the.taoofmac.com/space/blog/2005/10/11/2359
"""

import sys
from xml.dom import minidom, Node

class GPXParser:
    """ Used for reading out track data from the GPX files that can be
    exported by many hand-held GPS units. """
    def __init__(self, filename):
        self.tracks = {}
        try:
            doc = minidom.parse(filename)
            doc.normalize()
        except:
            sys.stderr.write('Error occurred parsing {0}\n'.format(filename))
        gpx = doc.documentElement
        for node in gpx.getElementsByTagName('trk'):
            self.ParseTrack(node)

    def ParseTrack(self, trk):
        name = trk.getElementsByTagName('name')[0].firstChild.data
        if not name in self.tracks:
            self.tracks[name] = {}
        for trkseg in trk.getElementsByTagName('trkseg'):
            for trkpt in trkseg.getElementsByTagName('trkpt'):
                lat = float(trkpt.getAttribute('lat'))
                lon = float(trkpt.getAttribute('lon'))
                ele = float(trkpt.getElementsByTagName('ele')[0].firstChild.data)
                time = trkpt.getElementsByTagName('time')[0].firstChild.data
                self.tracks[name][time]={'lat':lat,'lon':lon,'ele':ele}

    def List(self):
        """ List all track names within GPXParser. """
        return [str(i) for i in self.tracks.keys()]

    def GetTrack(self, name):
        """ Return a list of triples containing track coordinates, sorted by
        time. """
        times = self.tracks[name].keys()
        times.sort()
        points = [self.tracks[name][time] for time in times]
        return [(point['lon'],point['lat'],point['ele']) for point in points]

    def GetPolyline(self):
        pass

    def SplitNaturalBreaks(self):
        pass
