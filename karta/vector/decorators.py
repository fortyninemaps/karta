
def cache_decorator(key):
    """ Returns a method decorator that stores the result of a function
    invocation under `self._cache[key]`

    Does not update cache when *crs* is specified. """
    def wrapping_func(f):
        def replacement_func(self, *args, **kwargs):
            if (key in self._cache and
                    (args == self._cache[key]["args"]) and
                    (kwargs == self._cache[key]["kwargs"])):
                ret = self._cache[key]["output"]
            else:
                ret = f(self, *args, **kwargs)
                self._cache[key] = dict(output=ret, args=args, kwargs=kwargs)
            return ret
        return replacement_func
    return wrapping_func

