#!/usr/bin/env python

# Copyright (C) 2014-2017 Peter Batzel and Jason Sydes
#
# This file is part of Prost!.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the License with this program.
#
# Written by Jason Sydes, Peter Batzel.
# Advised by Thomas Desvignes.

# This file allows you to run prost (for development purposes) like so:
#   git clone ...prost
#   cd prost
#   python -m prost

# Python 3 imports
from __future__ import absolute_import
from __future__ import division

# Prost imports
from prost import prost

def main():
    prost.main()

if __name__ == "__main__":
    main()

# vim: softtabstop=4:shiftwidth=4:expandtab
