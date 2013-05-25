from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("raster.fill_sinks", ["raster/fill_sinks.pyx"]),
                   Extension("raster.cfuncs", ["raster/cfuncs.pyx"]),
                   Extension("vector._cvectorgeo", ["vector/_cvectorgeo.pyx"])]
)
