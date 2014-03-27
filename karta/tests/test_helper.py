""" Helper functions for testing """

import os
import inspect
import hashlib

TESTDIR = os.path.dirname(os.path.abspath(
                inspect.getfile(inspect.currentframe())))
TESTDATA = os.path.join(TESTDIR, "reference_data")

def md5sum(fnm, block=32768):
    """ Generate an MD5 hash from a file given by filename *fnm*. """
    md5 = hashlib.md5()
    with open(fnm, 'r') as f:
        while True:
            data = f.read(block)
            if not data:
                break
            md5.update(data)
    return md5.digest()
