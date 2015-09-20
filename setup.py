from ez_setup import use_setuptools
use_setuptools()
from os.path import exists
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext

VERSION = "0.5.0b1"

class build_ext(_build_ext):

    # solution taken from http://stackoverflow.com/questions/19919905/ \
    #   how-to-bootstrap-numpy-installation-in-setup-py
    def finalize_options(self):
        """ Unsets __NUMPY_SETUP__ so that the get_include function can be used
        """
        _build_ext.finalize_options(self)

        # Hack: unset __NUMPY_SETUP__ so that it can be imported
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())

        # Add include_dirs to each extension so that Cython will find them
        for extension in self.extensions:
            extension.include_dirs = self.include_dirs

        # Attempt to compile with Cython - or fallback on local C sources
        try:
            from Cython.Build import cythonize
            for ext in self.extensions:
                sources = []
                for src in ext.sources:
                    if src.endswith(".pyx"):
                        sources.append(src)
                    else:
                        sources.append(src+".pyx")
                ext.sources = sources

            _needs_stub = [ext._needs_stub for ext in self.extensions]
            self.extensions = cythonize(self.extensions)
            for ns, ext in zip(_needs_stub, self.extensions):
                ext._needs_stub = ns

        except ImportError:
            print("Cython not imported")
            print("If C sources exist in the source tree, they will be used")

            for ext in self.extensions:
                ext.sources = list(map(lambda f: f+".c", ext.sources))

        # Check that sources exist for each extension, either from Cython for
        # from having the C-code sitting around. Also, add a _needs_stub
        # attribute to each extension for setuptools.
        for ext in self.extensions:
            if not all(exists(src) for src in ext.sources):
                raise Exception("Extension sources not found - Cython required for install")
        return

# File extension is added to sources at overloaded build_ext.run()
extensions = [Extension("karta.raster.crfuncs", ["karta/raster/crfuncs"]),
              Extension("karta.vector._cvectorgeo", ["karta/vector/_cvectorgeo"])]

setup(
    name = "karta",
    version = VERSION,
    setup_requires = ["numpy>=1.6"],
    install_requires = ["numpy>=1.6", "pyproj>=1.9", "pyshp>=1.2"],
    author = "Nat Wilson",
    author_email = "njwilson23@gmail.com",
    packages = ["karta", "karta.vector", "karta.raster"],
    url = "http://www.ironicmtn.com/karta.html",
    description = "Geospatial analysis in Python",
    long_description = """
Karta - tidy package for geospatial computation
===============================================

*Karta* is a simple and fast framework for spatial analysis in Python.

Components:

- Clean geographically-aware vector and gridded data types
- Integration with pyproj to support a wide range of coordinate systems and
  transformations
- A selection of geographical analysis methods including geodetic length and
  area calculations, intersections, convex hulls, raster sampling and profiling,
  and grid warping
- IO for several common geographical formats, including GeoJSON, shapefiles
  (through pyshp), ESRI ASCII, and GeoTiff (through GDAL). Vector geometries
  implement ``__geo_interface__``.

*Karta* works with Python 2.7 and Python 3.3+.

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
    classifiers = ["Programming Language :: Python :: 2",
                   "Programming Language :: Python :: 2.7",
                   "Programming Language :: Python :: 3",
                   "Programming Language :: Python :: 3.2",
                   "Programming Language :: Python :: 3.3",
                   "Programming Language :: Python :: 3.4",
                   "Development Status :: 4 - Beta",
                   "Topic :: Scientific/Engineering"],
    license = "MIT License",
    ext_modules = extensions,
    cmdclass = {"build_ext": build_ext},
)
