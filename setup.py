from ez_setup import use_setuptools
use_setuptools()
from os.path import exists
from setuptools import setup, Extension
import numpy

VERSION = "0.4.0.1"

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

extensions = [Extension("karta.raster.cfill_sinks",
                        ["karta/raster/cfill_sinks"+ext],
                        include_dirs=include_dirs),
              Extension("karta.raster.crfuncs",
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

*Karta* is a Python/Python3 package for geospatial data structures and
computation. It provides simple and clean vector and raster data types, a
selection of geographical analysis methods, and the ability to read and write
several formats, including GeoJSON, shapefiles, and ESRI ASCII.

The goal of *Karta* is to expose a simple and fast framework for spatial
analysis.

DOCUMENTATION
-------------

See the `online manual`_.

.. _a link: http://www.ironicmtn.com/karta/doc/manual/karta-manual.html

The manual can also be built offline using Sphinx by running ``make`` from the
``doc/`` directory. The documentation is built from source code docstrings and
information in the `Wiki`_.

.. _a link: https://github.com/njwilson23/karta/wiki/Tutorial
""",
    download_url = "https://github.com/njwilson23/karta/archive/master.zip",
    classifiers = ["Programming Language :: Python",
                   "Programming Language :: Python :: 3",
                   "Development Status :: 4 - Beta",
                   "Topic :: Scientific/Engineering"],
    license = "MIT License",
    ext_modules = extensions
)

