#!/usr/bin/env bash

# setup autodoc
echo ; echo '--> Setting up autodoc...'
sphinx-apidoc -f -o docs/source prost

# Build docs
echo ; echo '--> Building docs...'
(cd docs; make html)

# Build it
echo ; echo '--> Building package...'
python setup.py build
