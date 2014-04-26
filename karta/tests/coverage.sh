#! /bin/bash
coverage run --source karta --omit *shapefile.py test_runner.py
coverage html
