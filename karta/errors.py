
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

class NoIntersection(Exception):
    """ Exception to raise when segments have no intersection. """
    def __init__(self, message=''):
        self.message = message

class GridError(Exception):
    def __init__(self, message=''):
        self.message = message
    def __str__(self):
        return self.message

class GridIOError(GridError):
    def __init__(self, message=''):
        self.message = message

class NonEquivalentGridError(GridError):
    def __init__(self, A, B, message=''):
        if len(message) == 0:
            self.message = ("{0} and {1} do not share equivalent "
                            "grid layouts".format(A, B))
        else:
            self.message = message

class MissingDependencyError(Exception):
    def __init__(self, message=''):
        self.message = message
    def __str__(self):
        return self.message

