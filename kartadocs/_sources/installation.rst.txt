Installation
============

The easiest way to install is to use ``pip``, although on some
systems ``conda`` may be easier for installing certain dependencies.
Installation requires ``setuptools``::

    pip install -U setuptools
    pip install karta

Building from source requires Cython and a C99-compliant compiler::

    pip install Cython

Clone the repository and build::

    git clone https://github.com/fortyninemaps/karta.git karta
    cd karta/
    python setup.py build

