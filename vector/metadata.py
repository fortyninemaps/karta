""" Metadata management """

# Metadata objects should
#   - maintain internal consistency
#   - be concatenable
#   - be iterable

import copy

class GeoMetadata(object):
    """ Class for handling collections of metadata. Data are organized similar
    to a Python dict, however different values must all have the same length
    and be of uniform type. Data are accessed by field (key). """

    _data = {}
    _fieldtypes = []

    def __init__(self, data):
        """ Create a collection of metadata from *data*, which may be a list
        with uniform type or a dictionary with equally-sized fields of uniform
        type. """

        if hasattr(data, 'keys') and hasattr(data.values, '__call__'):
            # Dictionary of attributes
            for k in data:
                dtype = type(data[k][0])
                if False in (isinstance(a, dtype) for a in data[k]):
                    raise GMetadataError("Data must have uniform type")

            n = len(data[k])
            if False in (len(data[k]) == n for k in data):
                raise GMetadataError("Data must have uniform lengths")

        else:
            # Single attribute
            if not hasattr(data, '__iter__'):
                data = [data]
            dtype = type(data[0])
            if False in (isinstance(a, dtype) for a in data):
                raise GMetadataError("Data must have uniform type")
            else:
                data = {'values': data}

        self._data = data
        self._fieldtypes = [type(data[k][0]) for k in data]

    def __add__(self, other):
        if isinstance(other, type(self)):
            res = copy.deepcopy(self)
        return res.extend(other)

    def __iter__(self):
        return self._data.__iter__()

    def __getitem__(self, idx):
        return tuple([self._data[k][idx] for k in self._data])

    def __delitem__(self, idx):
        for k in self._data:
            del self._data[k][idx]

    def getfield(self, name):
        """ Return all values from field *name*. """
        return self._data[name]

    def extend(self, other):
        """ Extend collection in-place. """
        if isinstance(other, type(self)):
            for i, k in enumerate(self._data):
                if type(self._fieldtypes[i]) == type(other._fieldtypes[i]):
                    self._data[k] += other._data[k]
                else:
                    raise GMetadataError("Cannot combine metadata instances "
                                         "with different type hierarchies")
        return self

    def keys(self):
        """ Return internal dictionary keys. """
        return self._data.keys()

    def values(self):
        """ Return internal dictionary values. """
        return self._data.values()


class GMetadataError(Exception):
    def __init__(self, message=''):
        self.message = message
    def __str__(self):
        return self.message


