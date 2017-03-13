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

# Python 3 imports
from __future__ import absolute_import
from __future__ import division

import sys

# For timing
import time

# For exception reporting/logging
import traceback

# Minimum python version
VERSION_PYTHON_MINIMUM =                                (2, 7)

def main_wrapper():
    """Wraps main() in some handy niceness."""

    # Get the start time.
    time_start = time.time()

    # Check python version first.
    if sys.version_info < VERSION_PYTHON_MINIMUM:
        print "Prost! requires Python {0} or later.  You are running {1}.".format(
            '.'.join([str(i) for i in VERSION_PYTHON_MINIMUM]),
            '.'.join([str(i) for i in sys.version_info[0:3]]))
        sys.exit(1)
    if sys.version_info >= (3,):
        print "Prost! doesn't support Python 3."
        sys.exit(1)

    # Now you can safely do Prost imports
    import prost.prost
    from prost.common import (ConfigurationException,
                              CannotContinueException,
                              fmt_time,
                              pmsg,
                              perr)

    try:
        prost.prost.main()
    except ConfigurationException as e:
        sys.stderr.write("\nConfiguration Error:\n{}\n\n".format(e))
        exit(1)
    except CannotContinueException as e:
        perr("\nError:\n{}\n\n".format(e))
    except KeyboardInterrupt as e:
        # We re-raise keyboard interrupts immediately.
        time_end = time.time()
        running_time = time_end - time_start
        pmsg("Total Prost running time: {}.".format(fmt_time(running_time)))
        raise e
    except Exception as e:
        perr("")
        perr(traceback.format_exc())

    time_end = time.time()
    running_time = time_end - time_start
    pmsg("Total Prost running time: {}.".format(fmt_time(running_time)))

# vim: softtabstop=4:shiftwidth=4:expandtab
