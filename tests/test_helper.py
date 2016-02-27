""" Defines some path variables for testing """

import os
import inspect

TESTDIR = os.path.dirname(os.path.abspath(
                inspect.getfile(inspect.currentframe())))
TESTDATA = os.path.join(TESTDIR, "data")
TMPDATA = os.path.join(TESTDIR, "tmp")
