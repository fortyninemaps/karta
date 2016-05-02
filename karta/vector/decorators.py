
def cache_decorator(key):
    """ Returns a method decorator that stores the result of a function
    invocation under `self._cache[key]`

    Does not update cache when *crs* is specified. """
    def wrapping_func(f):
        def replacement_func(self, *args, **kwargs):
            if key in self._cache:
                ret = self._cache[key]
            else:
                ret = f(self, *args, **kwargs)
                if "crs" not in kwargs:
                    self._cache[key] = (ret)
            return ret
        return replacement_func
    return wrapping_func

