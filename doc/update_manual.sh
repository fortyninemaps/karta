#! /bin/bash
#
# assumes that up to date RST sources are in build/html
# create these from the master branch, using
#
# ipython nbconvert tutorial.ipynb --to rst

cp build/html/* manual/ -ruv

