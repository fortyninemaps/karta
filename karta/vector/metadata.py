""" Metadata tables for vector data """

from collections import Sequence

class Metadata(Sequence):

    def __init__(self, data, fields=None, checktypes=False):
        """ Create a collection of metadata from *data*.

        *data* may be:
        
            - a list with uniform type
            - a dictionary with equally-sized fields of uniform type.
            - a scalar
            - another Metadata instance

        The *kwarg* fields is used for fast Metadata object creation, and
        requires *data* to be provided as a list of equal-sized tuples.

        If *checktypes* is True (default False), the values in *data* will be
        tested to ensure they are of constant type. In the usual case, this
        type-checking is foregone for speed.
        """
        if fields is None:
            if hasattr(data, "_fields"):        # Data is Metadata-like
                self._fields = data._fields
                self._data = data._data
            elif hasattr(data, "keys"):         # Data is dict-like
                self._fields = tuple(data.keys())
                fst = data[self._fields[0]]
                if hasattr(fst, "__iter__") and not isinstance(fst, str):   # Vector entries
                    zipvalues = zip(*[data[f] for f in self._fields])
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
            if hasattr(data[0], "__len__") and not isinstance(data[0], str):
                if len(data[0]) != len(fields):
                    raise ValueError("Length of data entries and fields don't match")
            self._fields = tuple(fields)
            self._data = data

        if checktypes:
            self._validatetypes()
        return

    def _validatetypes(self):
        for i,t in enumerate(self.types):
            if not all(isinstance(t, d[i]) for d in self._data):
                raise TypeError("data contains item ({0}) not of type {1}" \
                                .format(self._fields[i], t))
        return

    def __repr__(self):
        return "D[" + ", ".join(self._fields) + "]"

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

    def __len__(self):
        return len(self._data)

    def __contains__(self, field):
        return field in self._fields

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
        """ Return a dictionary for a single entry """
        d = self._data[i]
        if isinstance(i, slice):
            r = {}
            for j in range(len(self._fields)):
                r[self._fields[j]] = [d_[j] for d_ in d]
            return r
        else:
            return dict((self._fields[j], d[j]) for j in range(len(self._fields)))

    def getfield(self, field):
        """ Return list of data corresponding to *field* """
        if field in self._fields:
            i = self._fields.index(field)
            return [d[i] for d in self._data]
        else:
            raise KeyError("'{0}' not a field".format(field))

    def setfield(self, field, values):
        """ Modify or add a field with *values* """
        if len(values) != self.__len__():
            raise ValueError("mismatch between metadata length and field values")
        if field in self._fields:
            idx = self._fields.index(field)
            self._data = [tuplemut(row, v, idx) for row, v in zip(self._data, values)]
        else:
            idx = len(self._fields)
            self._fields = tupleinsert(self._fields, field, idx)
            self._data = [tupleinsert(row, v, idx) for row, v in zip(self._data, values)]
        return

    def extend(self, other):
        """ Extend Metadata from another Metadata instance. If *other* has
        field f in *self*, it is copied. Otherwise, the None is appended. """
        # TODO: appended None value should be approariate to the type of field f.
        idxs_other = []
        for f in self._fields:
            if f in other._fields:
                idxs_other.append(other._fields.index(f))
            else:
                idxs_other.append(None)

        for i in range(len(other)):
            self._data.append(tuple([None if idxs_other[j] is None
                                          else other._data[i][idxs_other[j]]
                                          for j in range(len(self._fields))]))

class Indexer(object):

    def __init__(self, md):
        if md is None:
            raise KeyError("cannot index data-less geometry")
        else:
            self.md = md

    def __getitem__(self, key):
        if isinstance(key, str):
            return self.md.getfield(key)
        elif isinstance(key, int):
            return self.md.get(key)
        else:
            raise KeyError("invalid key type: {0}".format(type(key)))

def tuplemut(tpl, val, idx):
    """ Return a tuple with *idx* changed to *val* """
    lst = list(tpl)
    lst[idx] = val
    return tuple(lst)

def tupleinsert(tpl, val, idx):
    """ Return a tuple with *val* inserted into position *idx* """
    lst = list(tpl[:idx]) + [val] + list(tpl[idx:])
    return tuple(lst)

