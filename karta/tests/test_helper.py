""" Helper functions for testing """

import os
import inspect
import hashlib
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

TESTDIR = os.path.dirname(os.path.abspath(
                inspect.getfile(inspect.currentframe())))
TESTDATA = os.path.join(TESTDIR, "reference_data")

def md5sum(s, block=32768):
    """ Generate an MD5 hash from a file, string buffer, or filename given by
    *s*. """
    md5 = hashlib.md5()
    try:
        if not hasattr(s, 'read'):
            s = open(s, 'r', encoding='utf-8')
        while True:
            data = s.read(block).encode('utf-8')
            if not data:
                break
            md5.update(data)
    finally:
        s.close()
    return md5.digest()

def tobuffer(s, pos=0):
    """ Place string *s* in a StringIO instance, set to position *pos*. """
    sio = StringIO(s)
    sio.seek(pos)
    return sio
