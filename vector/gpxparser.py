"""
GPX Parser

Based on code by Rui Carmo, available at
http://the.taoofmac.com/space/blog/2005/10/11/2359
"""

import sys
import numpy as np
from xml.dom import minidom, Node

try:
    from pyproj import Proj
except ImportError:
    # No projection capabilities
    pass

class GPXParser:
    """ Used for reading out track data from the GPX files that can be
    exported by many hand-held GPS units. """
    def __init__(self, filename):
        self.filename = filename
        self.tracks = {}
        self.has_projector = False
        try:
            doc = minidom.parse(filename)
            doc.normalize()
        except:
            sys.stderr.write('Error occurred parsing {0}\n'.format(filename))
        gpx = doc.documentElement
        for node in gpx.getElementsByTagName('trk'):
            self.parsetrack(node)

    def project(self, **kwargs):
        """ Wrapper for libproj.4 Proj class. """
        try:
            self.Projector = Proj(**kwargs)
            self.has_projector = True
        except NameError:
            raise ImportError("no pyproj found")
            return
        self.tracks_proj = {}
        for name in self.list():
            for t in self.tracks[name].keys():
                lon = self.track[name][t]['lon']
                lat = self.track[name][t]['lat']
                ele = self.track[name][t]['ele']
                x, y = self.Projector(lon, lat)
                self.tracks_proj[name][t] = {'x':x, 'y':y, 'ele':ele}

    def parsetrack(self, trk):
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

    def getnames(self):
        """ List all track names within GPXParser. """
        return [str(i) for i in self.tracks.keys()]

    def gettrack(self, name):
        """ Return a list of triples containing track coordinates
        sorted by time for track with name. """
        times = self.tracks[name].keys()
        times.sort()
        points = [self.tracks[name][t] for t in times]
        return [(point['lon'],point['lat'],point['ele']) for point in points]

    def getxyz(self, name):
        """ Return a list of triples containing track projected
        coordinates, sorted by time, for track with name. """
        times = self.tracks_proj[name].keys()
        times.sort()
        points = [self.tracks_proj[name][t] for t in times]
        return [(point['lon'],point['lat'],point['ele']) for point in points]

    def gettimes(self, name):
        """ Return a list of times for the track with name. """
        times = self.tracks[name].keys()
        times.sort()
        return times

    def getpolyline(self):
        """ Return a guppy.Polyline containing all tracks. """
        pass

    def filter_by_displacement(self, name, threshold, proj_kwargs=None):
        """ Return a list of lists, where each sublist contains consecutive GPS
        points between which the displacement is less than threshold. """
        # Make sure coordinates are projected
        if not self.has_projector:
            if proj_kwargs is not None:
                self.project(**proj_kwargs)
            else:
                raise ProjectionError("Need to define proj_kwargs dictionary")
                return
        xyz = self.gettrack_proj(name)
        X = np.array([i[0] for i in xyz])
        Y = np.array([i[1] for i in xyz])
        Z = np.array([i[2] for i in xyz])
        DIST = np.sqrt((X[:-1]-X[1:])**2 + (Y[:-1]-Y[1:])**2)
        #[INCOMPLETE - NEED TO FILTER BASED ON DIST]

    def split_natural_breaks(self, name):
        """ Split a track into linear natural segments, with nodes where
        direction changes. """
        pass


class ProjectionError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return self.value
