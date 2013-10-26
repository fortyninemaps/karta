from distutils.core import setup
from distutils.extension import Extension

try:
    from Cython.Build import cythonize
    USE_CYTHON = True
except ImportError:
    USE_CYTHON = False

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
    ext_modules = extensions
)

