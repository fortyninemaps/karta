Installation
============

The easiest way to install is to use ``pip``. Installation requires a
version of ``setuptools>=0.7.0``.

::

    pip install -U setuptools

To install the latest release from PyPI, run

::

    pip install karta

To build from source,

::

    git clone https://github.com/fortyninemaps/karta.git
    pip install -r karta/requirements.txt
    pip install karta/

Dependencies
------------

Required
~~~~~~~~

-  Python 2.7+ or Python 3.3+
-  numpy
-  pyshp
-  pyproj (for geodetic calculations)

Recommended
~~~~~~~~~~~

-  cython
-  gdal (for geotiff I/O)
-  scipy

When installing from PyPI, Cython-compiled C source code is provided and
will be automatically compiled to improve performance if a suitable C
compiler is available.


