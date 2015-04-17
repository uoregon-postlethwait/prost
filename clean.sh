#!/usr/bin/env bash

# clean docs
echo ; echo '--> Cleaning docs...'
(cd doc && make clean)

# setup version
echo ; echo '--> Cleaning build...'
python setup.py --verbose clean
