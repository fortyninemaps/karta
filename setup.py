from distutils.core import setup
from distutils.extension import Extension

try:
    import numpy
    from Cython.Build import cythonize
    include_dirs = [numpy.get_include()]
    USE_CYTHON = True
except ImportError:
    USE_CYTHON = False
    include_dirs = []
    print("Warning: not building accellerated modules")

ext = ".pyx" if USE_CYTHON else ".c"

extensions = [Extension("karta.raster.cfill_sinks", ["karta/raster/cfill_sinks"+ext]),
              Extension("karta.raster.crfuncs", ["karta/raster/crfuncs"+ext]),
              Extension("karta.vector._cvectorgeo", ["karta/vector/_cvectorgeo"+ext])]

if USE_CYTHON:
    extensions = cythonize(extensions)

setup(
    name = "karta",
    version = "0.1",
    author = "Nat Wilson",
    packages = ["karta", "karta.vector", "karta.raster"],
    include_dirs = include_dirs,
    ext_modules = extensions
)

