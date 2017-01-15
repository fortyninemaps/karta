
class GeometryError(Exception):
    """ Indicates invalid geometrical operations """
    def __init__(self, message=''):
        self.message = message
    def __str__(self):
        return self.message

class CRSError(Exception):
    """ Indicates invalid CRS parameters or operations """
    def __init__(self, message=''):
        self.message = message

class NoIntersection(Exception):
    """ Indicates that no intersection exists between segments """
    def __init__(self, message=''):
        self.message = message

class GridError(Exception):
    """ Indicates invalid Grid initialization, referencing, or transformation """
    def __init__(self, message=''):
        self.message = message
    def __str__(self):
        return self.message

