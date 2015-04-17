#!/usr/bin/env bash

# Build
./build.sh

# Link development egg
echo ; echo '--> Linking development egg...'
python setup.py develop
