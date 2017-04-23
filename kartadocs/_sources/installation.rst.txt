Installation
============

The easiest way to install in production is to use ``pip``, although on some
systems ``conda`` may be easier for installing certain dependencies.
Installation requires ``setuptools``::

    pip install -U setuptools

To install the latest release from PyPI, run::

    pip install karta


Building requires Cython and a C99-compliant compiler::

    pip install Cython

Clone the repository and build::

    git clone https://github.com/fortyninemaps/karta.git karta
    cd karta/
    python setup.py build

