
class GeometryError(Exception):
    """ Base class for geometry module errors. """
    def __init__(self, message=''):
        self.message = message
    def __str__(self):
        return self.message

class GInitError(GeometryError):
    """ Exception to raise when a geometry object fails to initialize. """
    def __init__(self, message=''):
        self.message = message
    def __str__(self):
        return self.message

class GUnitError(GeometryError):
    """ Exception to raise there is a projected unit problem. """
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

class MissingDependencyError(Exception):
    def __init__(self, message=''):
        self.message = message
    def __str__(self):
        return self.message

