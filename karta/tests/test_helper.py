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
    while True:
        data = s.read(block)
        data = data.encode('utf-8')
        print(data)
        if not data:
            break
        md5.update(data)
    return md5.digest()

def md5sum_file(fnm, block=32768):
    with open(fnm, 'r') as f:
        return md5sum(f, block=block)

def tobuffer(s, pos=0):
    """ Place string *s* in a StringIO instance, set to position *pos*. """
    sio = StringIO(s)
    sio.seek(pos)
    return sio
