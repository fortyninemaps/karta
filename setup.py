from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    name = "karta",
    version = "0.1",
    author = "Nat Wilson",
    #package_dir = {"karta": "src"},
    packages = ["karta", "karta.vector", "karta.raster"],
    cmdclass = {"build_ext": build_ext},
    ext_modules = [Extension("raster.cfill_sinks", ["karta/raster/cfill_sinks.pyx"]),
                   Extension("raster.crfuncs", ["karta/raster/crfuncs.pyx"]),
                   Extension("vector._cvectorgeo", ["karta/vector/_cvectorgeo.pyx"])],
)
