import sys
from ez_setup import use_setuptools
use_setuptools()
from os.path import exists
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext

sys.path.append("karta")
from version import __version__ as VERSION

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

        # Attempt to compile with Cython. Fallback on local C sources.
        #
        # Check that sources exist for each extension, either from Cython for
        # from having the C-code sitting around.
        #
        # Add a _needs_stub attribute to each extension for setuptools.
        try:
            from Cython.Build import cythonize
            _needs_stub = [ext._needs_stub for ext in self.extensions]
            self.extensions = cythonize(self.extensions)
            for ns, ext in zip(_needs_stub, self.extensions):
                ext._needs_stub = ns

        except ImportError:
            print("Cython not imported")
            print("If C sources exist in the build directory, they will be used")

            for ext in self.extensions:
                ext.sources = [f.replace(".pyx", ".c") for f in ext.sources]

        for ext in self.extensions:
            for src in ext.sources:
                if not exists(src):
                    raise Exception("Extension source code not found ({0})\n"
                                    "Cython required for install".format(src))
        return

# File extension is added to sources at overloaded build_ext.run()
extensions = [Extension("karta.raster.crfuncs", ["karta/raster/crfuncs.pyx"]),
              Extension("karta.vector.coordstring", ["karta/vector/coordstring.pyx"]),
              Extension("karta.vector.vectorgeo", ["karta/vector/vectorgeo.pyx"]),
              Extension("karta.vector.dateline", ["karta/vector/dateline.pyx"]),
              Extension("karta.vector.intersection", ["karta/vector/intersection.pyx"]),
              Extension("karta.vector.convexhull", ["karta/vector/convexhull.pyx"]),
              Extension("karta.vector.contains", ["karta/vector/contains.pyx"]),
              Extension("karta.vector.quadtree", ["karta/vector/quadtree.pyx"],
                        extra_compile_args=["-std=c99"]),
              Extension("karta.vector.rtree", ["karta/vector/rtree.pyx"],
                        extra_compile_args=["-std=c99"])]

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
    long_description = """*Karta* is a package for spatial analysis in Python. It streamlines data
processing by providing efficient generic geographical types for vector
and raster sources as well as a selection of analysis functions.

For example, read or create vector geometries:

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

Load and manipulate raster data:

.. code:: python

    grid = read_gtiff("landsat_scene.tif")  # Leverages GDAL
    grid.profile(line)              # Collect data along a line
    grid.resample(500.0, 500.0)     # Return a grid resampled at a new resolution

Documentation
-------------

See the `online
manual <http://www.fortyninemaps.com/kartadocs/introduction.html>`__,
the
`tutorial <http://www.fortyninemaps.com/kartadocs/_static/tutorial.html>`__,
or read the `API
documentation <http://www.fortyninemaps.com/kartadocs/reference.html>`__.

The manual can be built offline using
`Sphinx <http://sphinx-doc.org/>`__. Building the documentation requires
`numpydoc <https://github.com/numpy/numpydoc>`__.

Package Overview
----------------

-  **karta.crs**: framework for coordinate reference systems and
   geodetic calculations

-  **karta.vector.geometry**: geometry classes ``Point``,
   ``Multipoint``, ``Line``, and ``Polygon`` with associated methods
   such as length, area, intersections, membership testing, convex
   hulls, and affine transformations

-  **karta.raster.grid**: ``Grid`` classes including ``RegularGrid``
   class (supporting CRS-aware clipping, sampling, profiling along
   vector tracks), and experimental ``WarpedGrid``

-  **tests**: unit tests, to be run with ``python tests/runtests.py``

Formats
-------

*Karta* provides a basic working interface to several of common file
formats. Currently implemented are:

-  vector

   -  GeoJSON (r,w)
   -  ESRI Shapefiles (via GDAL) (r,w)
   -  GPS eXchange (GPX) (r,w)

-  raster

   -  GeoTiff (via GDAL) (r,w)
   -  ESRI ASCII Grid (r,w)

*Karta* implements the Python ```__geo_interface__``
attribute <https://gist.github.com/sgillies/2217756>`__ for vector
geometries. This permits data to be exchanged between *Karta* and
external modules that also implement ``__geo_interface__`` (e.g.
`shapely <https://github.com/Toblerity/Shapely>`__,
`fastkml <https://fastkml.readthedocs.org/en/latest/>`__).

Installation
------------

The easiest way to install in production is via ``pip``. Installation
requires a recent version of ``setuptools``:

::

    pip install -U setuptools

Then, to install the latest release from PyPI:

::

    pip install karta

Building from source
~~~~~~~~~~~~~~~~~~~~

Building from source requires Cython:

::

    pip install Cython

Then, clone the repository and install:

::

    git clone https://github.com/fortyninemaps/karta.git karta
    pip install -r karta/requirements.txt
    pip install karta/

Dependencies
------------

-  numpy >= >1.7
-  gdal >= 1.10
-  pyproj >= 1.9
-  blosc >= 1.2
-  C99-compliant compiler

*Karta* supports Python 2.7 and Python 3.4+.
""",
    classifiers = ["Programming Language :: Python :: 2",
                   "Programming Language :: Python :: 2.7",
                   "Programming Language :: Python :: 3",
                   "Programming Language :: Python :: 3.4",
                   "Programming Language :: Python :: 3.5",
                   "Topic :: Scientific/Engineering",
                   "Topic :: Scientific/Engineering :: GIS",
                   "License :: OSI Approved :: MIT License"],
    license = "MIT License",
    ext_modules = extensions,
    cmdclass = {"build_ext": build_ext},
)
