""" Helper functions for testing """

import hashlib

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
