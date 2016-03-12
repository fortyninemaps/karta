Installation
============

The easiest way to install in production is to use ``pip``, although on some
systems ``conda`` may be easier for installing certain dependencies.
Installation requires ``setuptools``.

::

    pip install -U setuptools

To install the latest release from PyPI, run

::

    pip install karta

Building from source uses *Cython*. Clone the repository and install:

::

    git clone https://github.com/fortyninemaps/karta.git karta
    pip install -r karta/requirements.txt
    pip install karta/


Dependencies
~~~~~~~~~~~~

- Python 2.7 or Python 3.3+
- numpy
- pyproj
- gdal
- C-compiler

