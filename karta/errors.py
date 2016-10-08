
class GeometryError(Exception):
    """ Base class for geometry module errors. """
    def __init__(self, message=''):
        self.message = message
    def __str__(self):
        return self.message

class GGeoError(GeometryError):
    """ Exception to raise when a geometry object attempts an invalid transform. """
    def __init__(self, message=''):
        self.message = message
    def __str__(self):
        return self.message

class CRSError(Exception):
    """ Exception to raise for invalid geodetic operations. """
    def __init__(self, message=''):
        self.message = message

class NoIntersection(Exception):
    """ Exception to raise when segments have no intersection. """
    def __init__(self, message=''):
        self.message = message

class GridError(Exception):
    def __init__(self, message=''):
        self.message = message
    def __str__(self):
        return self.message

