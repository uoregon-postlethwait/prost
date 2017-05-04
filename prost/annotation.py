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
# Written by Jason Sydes and Peter Batzel.
# Advised by Thomas Desvignes.

# Python 3 imports
from __future__ import absolute_import
from __future__ import division

# Prost imports
from prost.constants import *
from prost.common import (ExecutionException,
                          UnmappedAlignmentException,
                          pmsg,
                          perr,
                          SlotPickleMixin,
                          Progress)

# for abstract base classes:
from abc import ABCMeta, abstractmethod, abstractproperty

# Backported python 3.4 Enums, for annotations.  From package enum34.
from enum import Enum

# Pretty print the annotations
import pprint

###############
### Classes ###
###############


class Annotation(SlotPickleMixin):
    """Abstract Base Class for annotations."""

    __metaclass__ = ABCMeta
    __slots__ = ('name', 'kind')

    def __repr__(self):
        return pprint.pformat(self._key())

    def __str__(self):
        return self.name

    def _key(self):
        return (self.__class__.__name__, self.name, self.kind)

    def __eq__(self, other):
        return self._key() == other._key()

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self._key())


class MirbaseAnnotation(Annotation):
    """Abstract Base Class for miRBase annotations.

    Arguments:
        hit - An AlignmentExecutionHit object.
        species - String - The species under examination.
    """

    __metaclass__ = ABCMeta
    __slots__ = ()


class MirbaseMirAnnotation(MirbaseAnnotation):
    """A miRBase miR annotation."""
    __slots__ = ()

    class Kind(Enum):
        species = 0
        other_species = 1

    def __init__(self, hit, species):
        # Store the name (e.g. 'mmu-mir-205')
        self.name = hit.reference_sequence_name

        # Set the 'kind' based upon the species in the hit.
        # example reference sequence name in miRBase: 'mmu-mir-205'
        if hit.reference_sequence_name[:3] == species:
            self.kind = self.Kind.species
        else:
            self.kind = self.Kind.other_species


class MirbaseMirReverseAnnotation(MirbaseAnnotation):
    """A reverse miRBase miR annotation.

    In other words, derived from a 'reverse' annotation, where the query
    sequences are mature sequences from miRBase and the references is our
    compressed reads.
    """
    __slots__ = ()

    class Kind(Enum):
        species = 0
        other_species = 1

    def __init__(self, hit, species):
        # Store the name (e.g. 'mmu-mir-205')
        self.name = hit.query_sequence_name

        # Set the 'kind' based upon the species in the hit.
        # example reference sequence name in miRBase: 'mmu-mir-205'
        if hit.query_sequence_name[:3] == species:
            self.kind = self.Kind.species
        else:
            self.kind = self.Kind.other_species


class MirbaseHairpinAnnotation(MirbaseAnnotation):
    """A miRBase hairpin annotation."""
    __slots__ = ()

    class Kind(Enum):
        species = 0
        other_species = 1

    def __init__(self, hit, species):
        # Store the name (e.g. 'mmu-mir-205')
        self.name = hit.reference_sequence_name

        # Set the 'kind' based upon the species in the hit.
        # example reference sequence name in miRBase: 'mmu-mir-205'
        if hit.reference_sequence_name[:3] == species:
            self.kind = self.Kind.species
        else:
            self.kind = self.Kind.other_species


class BiomartOtherRNAAnnotation(Annotation):
    """..."""
    __slots__ = ('biotype',)

    # Because of peculiarities with the output, the only 'kind' of the
    # BiomartOtherRNAAnnotation is simply biomart_other_rna_annotation (in the
    # Mirbase annotations, the 'kind' corresponds to the column, and each
    # annotation fits inside one column; in the Biomart annotations, there is
    # only one kind, but each annotation is spread across two columns...)
    # (yah yah, kinda hacky...)
    class Kind(Enum):
        biomart_other_rna_annotation = 0

    def __init__(self, hit):

        # Store the kind (invariant) then name (e.g. Mir27b) then the biotype
        # of the annotation (e.g. miRNA, snoRNA, etc)
        self.kind = self.Kind.biomart_other_rna_annotation
        self.name = hit.reference_sequence_name.split('|')[0]
        self.biotype = hit.reference_sequence_name.split('|')[1]


# vim: softtabstop=4:shiftwidth=4:expandtab
