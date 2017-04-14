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
extensions = [
        Extension("karta.raster.crfuncs", ["karta/raster/crfuncs.pyx"]),

        Extension("karta.vector.coordstring", ["karta/vector/coordstring.pyx"]),

        Extension("karta.vector.vectorgeo", ["karta/vector/vectorgeo.pyx"],
                  extra_compile_args=["-std=c99"]),

        Extension("karta.vector.dateline", ["karta/vector/dateline.pyx"]),

        Extension("karta.vector.intersection", ["karta/vector/intersection.pyx"]),

        Extension("karta.vector.convexhull", ["karta/vector/convexhull.pyx"]),

        Extension("karta.vector.contains", ["karta/vector/contains.pyx"]),

        Extension("karta.vector.quadtree", ["karta/vector/quadtree.pyx"],
                        extra_compile_args=["-std=c99"]),

        Extension("karta.vector.rtree", ["karta/vector/rtree.pyx"],
                  extra_compile_args=["-std=c99"])
        ]

setup(
    name = "karta",
    version = VERSION,
    setup_requires = ["numpy>=1.10"],
    install_requires = ["numpy>=1.10",
                        "pyproj>=1.9",
                        "gdal>=1.10",
                        "blosc>=1.2.8",
                        "picogeojson>0.2.0"],
    author = "Nat Wilson",
    author_email = "natw@fortyninemaps.com",
    packages = ["karta", "karta.vector", "karta.raster"],
    url = "http://www.fortyninemaps.com/karta.html",
    description = "Geospatial analysis in Python",
    long_description = """
|Build Status| |Build status| |Coverage Status|

.. figure:: https://raw.githubusercontent.com/fortyninemaps/karta/gh-pages/images/karta_logo.png
   :alt: Karta

   Karta

*Karta* is a package for spatial analysis in Python. It simplifies
geospatial data processing by providing efficient generic classes for
vector and raster data sources, as well as a selection of analysis
functions.

Vector data types
-----------------

Data are represented as ``Point``, ``Line``, ``Polygon``,
``Multipoint``, ``Multiline``, and ``Multipolygon`` instances.

All data contain a ``.crs`` member encoding coordinate reference
information. All vector geometries possess a ``.properties`` dict
containing free-form metadata. *Multipart* geometries additionally
possess a ``.data`` member which is a simple typed table-like data
structure.

Geometries implement methods for computing distances, directions, and
spatial characteristics. *Multipart* geometries support fast spatial
indexing through quadtrees and r-trees.

GeoJSON and ESRI shapefile formats are supported for reading and
writing. Experimental support for GPX XML files is in the
``karta.vector.gpx`` submodule.

Vector geometries implement the Python ```__geo_interface__``
attribute <https://gist.github.com/sgillies/2217756>`__ for vector
geometries. This permits data to be exchanged between *Karta* and
external modules that also implement ``__geo_interface__`` (e.g.
`shapely <https://github.com/Toblerity/Shapely>`__,
`fastkml <https://fastkml.readthedocs.org/en/latest/>`__).

Raster data types
-----------------

The primary raster data type is the ``RegularGrid``, which represents
one or more 2-d arrays of pixels spaced via an affine transformation.
``RegularGrids`` are backed by one of several ``Band`` implementations,
with the default implementation using the
`blosc <http://www.blosc.org/>`__ compression library for efficient
in-memory storage. There is experimental support for disk-backed storage
via GDAL.

Grids may be queried, resampled, sliced, masked, and merged. Arbitrary
array-based functions may be mapped to raster data with
``RegularGrid.apply()``. Raster functions including slope, gradient, and
hillshade are in the ``karta.raster.misc`` submodule.

GeoTIFF images are the primary supported format, however ESRI ASCII
grids may also be used (with limitations due to the format).

Coordinate reference systems
----------------------------

Data in Karta is referenced to positions on earth via ``CRS`` objects
that implement projection and geodetic methods. Coordinate reference
systems may be either geographical or projected.

**Geographical CRS** objects return spatial relationships in terms of
the true computed distances and azimuths on a spherical or ellipsoidal
Earth.

**Projected CRS** objects (e.g. UTM, Polar Stereographic, and Web
Mercator) return spatial relationships in terms of a flat plane,
dependent on the projection.

Examples
--------

Read or create vector geometries:

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

Installation
------------

*Karta* currently supports Python 2.7 and Python 3.4+.

The easiest way to install is via ``pip``. Installation requires a recent
version of ``setuptools``.

::

    pip install -U setuptools
    pip install karta

Building from source
~~~~~~~~~~~~~~~~~~~~

Building from source requires Cython and a C99-compliant compiler:

::

    pip install Cython

Then, clone the repository and install:

::

    git clone https://github.com/fortyninemaps/karta.git karta
    cd karta/
    python setup.py build


Documentation
-------------

See the `online
manual <http://www.fortyninemaps.com/kartadocs/introduction.html>`__,
the
`tutorial <http://www.fortyninemaps.com/kartadocs/_static/tutorial.html>`__,
or read the `API
documentation <http://www.fortyninemaps.com/kartadocs/reference.html>`__.

Contributing
------------

Bug reports, feature requests, and pull requests are welcome.

Run unit tests with ``python tests/runtests.py``.

The manual is built using `Sphinx <http://sphinx-doc.org/>`__ and
requires `numpydoc <https://github.com/numpy/numpydoc>`__.

.. |Build Status| image:: https://travis-ci.org/fortyninemaps/karta.svg?branch=master
   :target: https://travis-ci.org/fortyninemaps/karta
.. |Build status| image:: https://ci.appveyor.com/api/projects/status/viiimwp5pu7ff2bp?svg=true
   :target: https://ci.appveyor.com/project/njwilson23/karta
.. |Coverage Status| image:: https://coveralls.io/repos/github/fortyninemaps/karta/badge.svg?branch=master
   :target: https://coveralls.io/github/fortyninemaps/karta?branch=master
""",
    classifiers = ["Programming Language :: Python :: 2",
                   "Programming Language :: Python :: 2.7",
                   "Programming Language :: Python :: 3",
                   "Programming Language :: Python :: 3.4",
                   "Programming Language :: Python :: 3.5",
                   "Programming Language :: Python :: 3.6",
                   "Topic :: Scientific/Engineering",
                   "Topic :: Scientific/Engineering :: GIS",
                   "License :: OSI Approved :: MIT License"],
    license = "MIT License",
    ext_modules = extensions,
    cmdclass = {"build_ext": build_ext},
)
