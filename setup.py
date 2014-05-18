import os
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
    print("Warning: Cython not imported.")
    print("If C sources exist in the source tree, they will be used.")

ext = ".pyx" if USE_CYTHON else ".c"

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
    if False in (os.path.exists(src) for src in extension.sources):
        # C-extension sources missing, so don't try to build them
        extensions = []
        print("Warning: Extension source not found.")
        print("Not building accelerated modules")
        break

setup(
    name = "karta",
    version = "0.1",
    author = "Nat Wilson",
    author_email = "njwilson23@gmail.com",
    packages = ["karta", "karta.vector", "karta.raster"],
    ext_modules = extensions
)

