from ez_setup import use_setuptools
use_setuptools()
from os.path import exists
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext

VERSION = "0.6.1"

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
            print("If C sources exist in the build directory, they will be used")

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
    install_requires = ["numpy>=1.6", "pyproj>=1.9", "gdal>=1.10", "blosc>=1.2.8"],
    author = "Nat Wilson",
    author_email = "njwilson23@gmail.com",
    packages = ["karta", "karta.vector", "karta.raster"],
    url = "http://www.fortyninemaps.com/karta.html",
    description = "Geospatial analysis in Python",
    long_description = """
*Karta* streamlines processing raster and vector geographical data in Python.

Create vector geometries:

.. code:: python

    point = Point((-130.0, 52.0), crs=LonLatWGS84)

    line = read_geojson("linedata.json")

    polygon = Polygon([(-515005.78, -1301130.53),
                       (-579174.89, -1282271.94),
                       (-542977.83, -1221147.82),
                       (-437864.05, -1251641.55),
                       (-438160.72, -1252421.48),
                       (-437961.28, -1285314.00)],
                       crs=NSIDCNorth)

Perform simple queries:

.. code:: python

    point2 = Point((-25.0, 48.0), crs=LonLatWGS84)
    point.distance(point2)          # Distance in geographical units

    line.intersects(polygon)        # True or False

    ch = polygon.convex_hull()      # Returns a new polygon
    ch.to_shapefile("poly.shp")

Work with raster data:

.. code:: python

    grid = read_gtiff("landsat_scene.tif")  # Leverages GDAL

    grid.profile(line)              # Collect data along a line

    grid.resample(500.0, 500.0)     # Return a grid resampled at a new resolution

*Karta* works with Python 2 and 3. Suggestions, bug reports, test cases, and
pull requests are welcome.

DOCUMENTATION
-------------

See the `webpage <http://www.fortyninemaps.com/karta.html>`__ and the `manual
<http://www.fortyninemaps.com/kartadocs/karta-manual.html>`__.

The manual can also be built offline using Sphinx by running ``make`` from the
``doc/`` subdirectory. The documentation is built from source code docstrings
and the example IPython notebooks, which are also reproduced in the `Wiki
<https://github.com/fortyninemaps/karta/wiki/Tutorial>`__. Building the
documentation requires `Sphinx <http://sphinx-doc.org/>`__, `alabaster
<https://github.com/bitprophet/alabaster>`__ and `numpydoc
<https://github.com/numpy/numpydoc>`__.

DEPENDENCIES
------------

Required
~~~~~~~~

- numpy
- blosc
- GDAL
- pyproj

When installing from PyPI, C source code is provided. When building from
sources, Cython is required.
""",
    classifiers = ["Programming Language :: Python :: 2",
                   "Programming Language :: Python :: 2.7",
                   "Programming Language :: Python :: 3",
                   "Programming Language :: Python :: 3.3",
                   "Programming Language :: Python :: 3.4",
                   "Programming Language :: Python :: 3.5",
                   "Topic :: Scientific/Engineering",
                   "Topic :: Scientific/Engineering :: GIS",
                   "License :: OSI Approved :: MIT License"],
    license = "MIT License",
    ext_modules = extensions,
    cmdclass = {"build_ext": build_ext},
)
