# Copyright (C) 2014, 2015  Peter Batzel and Jason Sydes
#
# This file is part of Prost!.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the License with this program.
#
# Written by Jason Sydes and Peter Batzel.
# Advised by Thomas Desvignes.

# Python 3 imports
from __future__ import absolute_import
from __future__ import division

import sys

# For progress and timing
import time

# For logging
import logging


##################
### Exceptions ###
##################

class ConfigurationException(Exception):
    """Something went wrong with the configuration section."""
class ExecutionException(Exception):
    """Something went wrong with the execution of an external command."""
class ControlFlowException(Exception):
    """Something went wrong with control flow logic."""
class CannotContinueException(Exception):
    """Prost has encountered a situation from which it cannot continue."""
class PrerequisitesException(Exception):
    """Prost is missing prerequisites (e.g. an out of date BBMap version)."""
class SamExtendedCigarParsingException(Exception):
    """Something when wrong when parsing a SAM Extended CIGAR string."""
class SamMDTagParsingException(Exception):
    """Something when wrong when parsing a SAM MD Tag."""
class MirModificationCalculationException(Exception):
    """Problems with the calculation of miRNA modification percents.  This
    exception exists to help us catch unexpected situations."""
class ModificationThingEncounteredNAlignment(Exception):
    """While building the ModifcationThing, either the main or member genomic
    location had a "M" in it, meaning that it aligned to an 'N' in the
    reference genome.  We ignore these for the calculation and log their
    existance."""
class ArgumentTypeException(Exception):
    """Passed in the incorrect type of argument."""
class UnmappedAlignmentException(Exception):
    """Generally the aligner should be configured to not output unmapped reads.
    In case it is outputting unmapped reads, raise (and catch) this exception.
    Originally added because in one case BBMapSkimmer was configured properly
    with outputunmapped=f, but it was outputting a single weird read.

    Example:
        TACTCAGGATCTAGTTTTCCAACTTTGC 272 17 50130514 24 * * 0 0 * * AM:i:24 NH:i:2 YR:i:2090 """


#################
### Functions ###
#################

def all_slots(cls):
    """Return a list of all __slots__ in the class hierarchy.

    Used by SlotPickleMixin.
    """

    def _class_list(cls):
        """Return a list which includes cls and superclasses of cls."""
        classes = [cls]
        for c in cls.__bases__:
            if c not in (object, SlotPickleMixin):
                classes += _class_list(c)
        return classes

    slots = []
    for c in _class_list(cls):
        if hasattr(c, '__slots__'):
            slots += c.__slots__
    return slots

def pmsg(msg):
    """Prost Message.  Write a message to STDOUT and the prost.log.

    You should *never* use print(), sys.stderr.write() or sys.stdout.write()
    within Prost (except for the class Progress).  Instead, use this function.

    """
    logging.info(msg)
    print(msg)

def perr(msg):
    """Prost Error Message.  Same as pmsg(), but it's an error instead of info."""

    logging.error(msg)
    print(msg)

def fmt_time(seconds):
    """Format number of seconds into a human readable time string.

    Example:

        Example::

            3601.5 ->  "3601.5s (or 1h0m)"
            1000.5 ->  "1000.5s (or 16m40s)"
            59.5 ->  "59.5s"

    Arguments:
        seconds (float): A length of time given in seconds.

    """
    m, s = divmod(seconds, 60)
    h, m = divmod(int(m), 60)
    s = int(s)
    if h >= 1:
        return "{0:.1f}s (or {1}h{2}m)".format(int(seconds), h, m)
    elif m >= 1:
        return "{0:.1f}s (or {1}m{2}s)".format(seconds, m, s)
    else:
        return "{0:.1f}s".format(seconds)


###############
### Classes ###
###############

class SlotPickleMixin(object):
    """Mixin to allow pickling of objects whose class defines __slots__.

    Taken from
    http://code.activestate.com/recipes/578433-mixin-for-pickling-objects-with-__slots__/

    Note that this does not currently work with inheritance, because __slots__
    only looks at the current class.  It's possible this could be made to work
    with inheritance, but it's probably easier to simply avoid inheritance and
    instead use mixins to provide common functionality (see GenomicLocation and
    GenomicLocationWithDesignation for an example; the latter was previously
    a subclass of the former, but was rewritten as just described).
    """

    def __getstate__(self):
        return dict(
            (slot, getattr(self, slot))
            for slot in all_slots(self.__class__)
            if hasattr(self, slot)
        )

    def __setstate__(self, state):
        for slot, value in state.items():
            setattr(self, slot, value)


class Progress(object):
    """A simple progress meter printing to STDERR and prost.log.

    We print to STDERR just for convenience because it's not buffered.

    """
    def __init__(self, progress_str, step=None, known_total=None, **kwargs):
        self.start = time.time()
        self.count = 0
        self.total = 0
        self.progress_str = progress_str
        self.known_total = None if known_total is None else int(known_total)
        self.step = None if step is None else int(step)
        self.indent = kwargs.get('indent', 0) * ' '
        sys.stderr.write("{}{}...".format(self.indent, progress_str))
        if self.step:
            sys.stderr.write('\n')

    def progress(self):
        if self.step:
            self.count += 1
            self.total += 1
            if self.count == self.step:
                if self.known_total:
                    sys.stderr.write("\r{}    {}/{} so far...".format(
                        self.indent, self.total, self.known_total))
                else:
                    sys.stderr.write("\r{}    {} so far...".format(
                        self.indent, self.total))
                self.count = 0
        else:
            sys.stderr.write('.')

    def done(self):
        _end = time.time()
        _secs = _end - self.start

        if not self.step:
            msg_mid = "{}{}...".format(self.indent, self.progress_str)
            sys.stderr.write('\r')
            sys.stderr.write(msg_mid)
            msg_mid = ""
        else:
            if self.known_total:
                msg_mid = "{}    {}/{} so far...".format(
                    self.indent, self.total, self.known_total)
                sys.stderr.write('\r')
                sys.stderr.write(msg_mid)
            else:
                msg_mid = "{}    {} so far...".format(self.indent, self.total)
                sys.stderr.write('\r')
                sys.stderr.write(msg_mid)

        msg_end = " done. (elapsed time: {})".format(fmt_time(_secs))
        sys.stderr.write(msg_end + '\n')

        # Now write full message to the log file.

        msg_beg = ("{}{}...".format(self.indent, self.progress_str))
        if self.step:
            msg_beg += '\n'
        logging.info(msg_beg + msg_mid + msg_end)


# vim: softtabstop=4:shiftwidth=4:expandtab
