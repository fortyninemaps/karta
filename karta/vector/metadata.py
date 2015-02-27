""" Faster Metadata implementation

Based on the observation that constructing Metadata is very common, while
changing it is less so.
"""

from collections import Sequence

class Metadata(Sequence):

    def __init__(self, data, fields=None, checktypes=False):
        if fields is None:
            if isinstance(data, Metadata):      # Data is a Metadata instance
                self._data = data._data
                self._fields = data._fields
            elif hasattr(data, "keys"):         # Data is dict-like
                self._fields = tuple(data.keys())
                fst = data[self._fields[0]]
                if hasattr(fst, "__iter__") and not isinstance(fst, str):   # Vector entries
                    zipvalues = zip(*(data[f] for f in self._fields))
                    self._data = [tuple(item) for item in zipvalues]
                else:                           # Scalar entries
                    self._data = [(data[f],) for f in self._fields]
            else:                               # Data is list or scalar
                self._fields = ("value",)
                if hasattr(data, "__iter__") and not isinstance(data, str):
                    self._data = [(d,) for d in data]
                else:
                    self._data = [(data,)]
        else:
            if hasattr(data[0], "__len__") and not isinstance(data[0], str):    # TODO: avoid assertions
                assert(len(data[0]) == len(fields))
            self._fields = tuple(fields)
            self._data = data

        if checktypes:
            self._validatetypes()
        return

    def _validatetypes(self):
        for i,t in enumerate(self.types):
            if not all(isinstance(t, d[i]) for d in self._data):
                raise TypeError("data contains item ({0}) not of type {1}".format(self._fields[i], t))
        return

    def __eq__(self, other):
        return (self._data == other._data) and (self._fields == other._fields)

    def __neq__(self, other):
        return self != other

    def __getitem__(self, i):
        return self._data[i]

    def __setitem__(self, i, val):
        self._data[i] = val
        return

    def __delitem__(self, i):
        del self._data[i]

    def __len__(self, i):
        return len(self._data)

    @property
    def data(self):
        return self._data

    @property
    def fields(self):
        return self._fields

    @property
    def types(self):
        return tuple(type(a) for a in self._data[0])

    def get(self, i):
        d = self._data[i]
        if isinstance(i, slice):
            r = {}
            for j in range(len(self._fields)):
                r[self._fields[j]] = [d_[j] for d_ in d]
            return r
        else:
            return dict((self._fields[j], d[j]) for j in range(len(self._fields)))

    def getfield(self, name):
        if name in self._fields:
            i = self._fields.index(name)
            return [d[i] for d in self._data]
        else:
            raise KeyError("'{0}' not a field name".format(name))

