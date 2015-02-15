from ez_setup import use_setuptools
use_setuptools()
from os.path import exists
from setuptools import setup, Extension
import numpy

VERSION = "0.4.2.1"

try:
    from Cython.Build import cythonize
    USE_CYTHON = True
    include_dirs = [numpy.get_include()]
    ext = ".pyx"
except ImportError:
    USE_CYTHON = False
    include_dirs = []
    ext = ".c"
    print("Warning: Cython not imported")
    print("If C sources exist in the source tree, they will be used")

extensions = [Extension("karta.raster.crfuncs",
                        ["karta/raster/crfuncs"+ext],
                        include_dirs=include_dirs),
              Extension("karta.vector._cvectorgeo",
                        ["karta/vector/_cvectorgeo"+ext],
                        include_dirs=include_dirs)]

if USE_CYTHON:
    extensions = cythonize(extensions)

for extension in extensions:
    if not all(exists(src) for src in extension.sources):
        # C-extension sources missing, so don't try to build them
        extensions = []
        print("Warning: Extension source not found: {0}".format(extension.sources))
        print("Not building accelerated modules")
        break

setup(
    name = "karta",
    version = VERSION,
    install_requires = ["numpy>=1.6", "cython>=0.15", "pyproj>=1.9", "pyshp>=1.2"],
    author = "Nat Wilson",
    author_email = "njwilson23@gmail.com",
    packages = ["karta", "karta.vector", "karta.raster"],
    url = "http://www.ironicmtn.com/karta.html",
    description = "Geospatial analysis in Python",
    long_description = """
Karta - tidy package for geospatial computation
===============================================

*Karta* provides a simple and fast framework for spatial analysis in Python.

Components:

- Clean vector and raster data types that are geographical coordinate
  system-aware
- A selection of geographical analysis methods including geodetic length and
  area calculations, intersections, convex hulls, raster sampling, and grid
  warping
- IO for several common geographical formats, including GeoJSON, shapefiles,
  and ESRI ASCII

*Karta* works with Python 2.6-2.7 and Python 3.3+.

DOCUMENTATION
-------------

See the `online manual <http://www.ironicmtn.com/kartadocs/karta-manual.html>`_,
read the tutorial_, or search the `API documentation`_.

.. _tutorial: http://www.ironicmtn.com/kartadocs/tutorial.html
.. _API documentation: http://www.ironicmtn.com/kartadocs/reference.html

The manual can also be built offline with Sphinx by running ``make`` from the
``doc/`` directory. The documentation is built from source code docstrings and
information in the `Wiki <https://github.com/njwilson23/karta/wiki/Tutorial>`_.
""",
    download_url = "https://github.com/njwilson23/karta/archive/master.zip",
    classifiers = ["Programming Language :: Python",
                   "Programming Language :: Python :: 3",
                   "Development Status :: 4 - Beta",
                   "Topic :: Scientific/Engineering"],
    license = "MIT License",
    ext_modules = extensions
)

