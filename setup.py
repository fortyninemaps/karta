from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("raster.cfill_sinks", ["raster/cfill_sinks.pyx"]),
                   Extension("raster.crfuncs", ["raster/crfuncs.pyx"]),
                   Extension("vector._cvectorgeo", ["vector/_cvectorgeo.pyx"])]
)
