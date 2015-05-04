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

# Prost imports
from prost.constants import *
from prost.common import (ExecutionException,
                          UnmappedAlignmentException,
                          pmsg,
                          perr,
                          SlotPickleMixin,
                          Progress)
#from prost.db_caching import (DB_PROXY,
#                              SequenceDB,
#                              ShortSeqCache,
#                              GenomicLocationCache)

import sys
import os.path
import shlex
import subprocess

# for logging
import logging

# for abstract base classes:
from abc import ABCMeta, abstractmethod, abstractproperty

# Converting strings to classes
import importlib

# For nicer output (e.g. BBMap command line)
import textwrap

# For parsing cigar strings
import re, itertools

# For quick element grabbing during construction of alignment execution hits
import operator

# For debugging
import traceback

# For pretty printing
import pprint

#################
### Functions ###
#################

def str_to_cls(module_name, class_name):
    """Convert a string to a python class object.

    Arguments:
        module_name (str): The string representation of the module containing
            the class.
        class_name (str): The string representation of the class.

    Returns:
        type: A python class object.

    """
    # load the module, will raise ImportError if module cannot be loaded
    m = importlib.import_module(module_name)
    # get the class, will raise AttributeError if class cannot be found
    c = getattr(m, class_name)
    return c


###############
### Classes ###
###############

class Alignments(object):
    """A singleton representing the Sequence Alignments to be performed.

    Assumes that the 'Configuration' object passed in has been sanitized."""

    @property
    def genome_alignment(self):
        """Returns the only genome alignment, which is the first alignment we
        store. See Alignments.annotation_alignments."""
        return self._alignments[0]

    @property
    def annotation_alignments(self):
        """Returns the annotation alignments, which are every alignment except
        the first alignment.  See Alignments.genome_alignment."""
        return self._alignments[1:]

    def __init__(self, config):
        """Initialization.

        Arguments:
            config: The 'Configuration' singleton object.
        """
        self._alignments = []
        self._create_alignments(config)

    def align(self, conf):
        """Perform the sequence alignments.

        Arguments:
            conf: The Configuration singleton object.
        """
        progress = Progress("Alignments")

        skip = (conf.general.skip_genome_alignment or
                conf.general.skip_sequence_alignments)
        self.genome_alignment.align(skip)
        skip = (conf.general.skip_annotation_alignments or
                conf.general.skip_sequence_alignments)
        for annotation_alignment in self.annotation_alignments:
            annotation_alignment.align(skip)

        progress.done()

    def _create_alignments(self, conf):

        # Create the Genome Alignment
        self._create_alignment(conf, 'GenomeAlignment', GenomeAlignment)

        # Create any Annotation Alignments
        for anno_align_name in conf.annotation_alignment_names:
            self._create_alignment(conf, anno_align_name, AnnotationAlignment)

    def _create_alignment(self, conf, conf_alignment_section_name, align_cls):
        """Create an Alignment object.

        Arguments:
            conf: The Configuration singleton object.
            conf_alignment_section_name: The Configuration section name
                corresponding to this Alignment.
            align_cls: The class of the Alignment object to create (e.g.
                GenomeAlignment)
        """

        # Set the Annotation class, if any.
        conf_alignment = conf.get_section(conf_alignment_section_name)
        if align_cls == AnnotationAlignment:
            annotation_cls = str_to_cls('prost.prost', conf_alignment['type'])

        # Create the Alignment object.
        #     Examples:
        #     self._alignments.append(GenomeAlignment(BBMapAlignmentExecution,
        #         'GenomeAlignment'))
        #     self._alignments.append(AnnotationAlignment(BBMapAlignmentExecution,
        #         'AnnotationAlignment2'))
        if align_cls == GenomeAlignment:
            self._alignments.append(align_cls(conf,
                conf_alignment_section_name))
        elif align_cls == AnnotationAlignment:
            self._alignments.append(align_cls(conf,
                conf_alignment_section_name, annotation_cls))
        else:
            raise ControlFlowException, \
                    "ERR541: Shouldn't be possible to reach here."


class Alignment(object):
    """Abstract Base Class for a type of alignment.

    Currently there are two types of sequence alignments we perform:
        1. An alignment to the genome (GenomeAlignment).
        2. Several alignments to produce annotations (AnnotationAlignment).

    This class delegates the responsibility of actually executing the search to
    the AlignmentExecution class.  Hence an Alignment object "has_one"
    AlignmentExecution object.
    """

    __metaclass__ = ABCMeta

    def __init__(self, conf, conf_alignment_section_name, input_file):
        """Creates the contained AlignmentExecution object.

        Arguments:
            conf: The Configuration singleton object.
            conf_alignment_section_name: The Configuration section name
                corresponding to this Alignment.
        """

        # Based upon the search tool specified in configuration, create the
        # AlignmentExecution object to use for this Alignment.
        conf_alignment = conf.get_section(conf_alignment_section_name)
        tool = conf_alignment['tool'].lower()
        if tool == 'bbmap':
            align_exec_cls = BBMapAlignmentExecution
        else:
            raise ControlFlowException, \
                    "ERR323: Shouldn't be possible to reach here."

        self.alignment_execution = align_exec_cls(conf,
                conf_alignment_section_name, input_file)

    def align(self, skip):
        """Performs the sequence alignment.
        Arguments:
            skip: Boolean - Whether or not to skip this particular alignment.
        """
        self.alignment_execution.align(skip)

    def write_to_cache(self, short_seqs, db_proxy):
        """This writes out the results of a sequence alignment to a database
        cache.  For each non-previously-cached short sequence, it saves the
        following information:

        * the string representation of the short sequence itself
        * its designation
        * either:
            1. the alignments hits comprising the designation as genomic
               locations (note that hits with lower designations are excluded)
            2. the fact that no alignment hit was found
        * the (target) sequence database used

        Arguments:
            short_seqs: The ShortSeqs singleton
            db: The database handle
        """

        ### NOTE -  Currently broken! We stopped using GenomicLocation and instead are now just using AlignmentExecutionHits.
        # Get the sequence database first
        db = SequenceDB.get(SequenceDB.name == self.alignment_execution.seq_db)

        # Cache it
        progress = Progress(
                "Writing to database cache, this may take some time",
                100, len(short_seqs))
        with db_proxy.transaction():
            for short_seq in short_seqs.itervalues():
                if short_seq.cached:
                    continue
                short_seq_cached = ShortSeqCache.create(db=db,
                        sequence=short_seq.seq_str,
                        designation_integer=short_seq.designation_integer)
                for gen_loc in short_seq.genomic_locations:
                    # Store the full designation if it's a 'three' or more.
                    if isinstance(gen_loc, GenomicLocationWithDesignation):
                        # # TODO: But don't cache no_hit sequences???
                        # if gen_loc.designation == DESIGNATION_NO_HIT:
                        #     continue
                        GenomicLocationCache.create(short_seq=short_seq_cached,
                                lg=gen_loc.lg, start=gen_loc.start,
                                end=gen_loc.end,
                                designation=gen_loc.designation)
                    else:
                        GenomicLocationCache.create(short_seq=short_seq_cached,
                                lg=gen_loc.lg, start=gen_loc.start,
                                end=gen_loc.end)
                short_seq.cached = True
                progress.progress()
        progress.done()
        # TODO: This may be very slow performance wise.  Could switch to
        # pickled gen_locs.  We'll see. Update: it does look pretty slow.
        # Alternative is to write extended inserts (would have to write them
        # specific to each database driver, as mysql, postgres, and presumably
        # sqlite do these all differently.... (e.g. mysql has extended inserts
        # while postgres has COPY).


class GenomeAlignment(Alignment):
    """ ... """
    def __init__(self, conf, conf_alignment_section_name):
        """

        Arguments:
            conf: The Configuration singleton object.
            conf_alignment_section_name: The Configuration section name
                corresponding to this Alignment.
        """

        input_file = conf.general.input_file_for_alignments
        super(GenomeAlignment, self).__init__(conf,
                conf_alignment_section_name, input_file)


class AnnotationAlignment(Alignment):
    """ ... """

    def __init__(self, conf, conf_alignment_section_name, annotation_cls):
        """Same as Alignment, except for one additional argument.

        Arguments:
            conf: The Configuration singleton object.
            conf_alignment_section_name: The Configuration section name
                corresponding to this Alignment.
            annotation_cls: The Annotation class used for annotating.
        """

        self.annotation_cls = annotation_cls
        input_file = conf.general.input_file_for_alignments
        super(AnnotationAlignment, self).__init__(conf,
                conf_alignment_section_name, input_file)


class AlignmentExecution(object):
    """Abstract Base Class for alignment executions (e.g. BBMap searches).

    Represents an actual execution of a sequence alignment.
    """

    __metaclass__ = ABCMeta

    @property
    def seq_db(self): return self._seq_db

    @property
    def parameters_file(self): return self._parameters_file

    @property
    def input_file(self): return self._input_file

    @property
    def results_file(self): return self._results_file

    @property
    def command_line(self):
        """Return an array of strings that, when joined into a single string,
        is the genomic search command that will be executed on the system.

        This is a read-only property. This is because the attribute is
        calculated from an overlay of defaults, the parameter file, and passed
        options.  See generate_command_line().
        """
        return self._cline

    def __init__(self, conf, conf_alignment_section_name, input_file):
        """
        Arguments:
            conf: The Configuration singleton object.
            conf_alignment_section_name: The name of the Alignment section
                corresponding to this Alignment.
        """

        # Get the Alignment configuration sections
        conf_alignment = conf.get_section(conf_alignment_section_name)

        # Set the results_file name.
        try:
            results_file = conf_alignment['results_file']
        except KeyError:
            # results_file wasn't specified. Set to default.
            results_file = "alignment_results_{}_{}".format(
                    conf_alignment['tool'].lower(), conf_alignment['name'])

        # Initialize
        self._cline = ""
        self._seq_db = conf_alignment['db']
        self._input_file = input_file
        self._results_file = results_file
        self._max_threads = conf.general.max_threads
        self._min_seq_length = conf.general.min_seq_length
        self.max_3p_mismatches = conf_alignment['max_3p_mismatches']
        self.max_non_3p_mismatches = conf_alignment['max_non_3p_mismatches']
        self.indelnt_penalty_multiplier = \
            conf_alignment.get('indelnt_penalty_multiplier',
                    DEFAULT_INDELNT_PENALTY_MULTIPLIER)
        self._allow_indels = conf_alignment['allow_indels']
        self.generate_command_line()

    @abstractmethod
    def generate_command_line(self):
        """Generate the command line of the search command.

        Arguments:
            override_commands: A dict, containing key/val command line args.
                key: command parameter (e.g. '-task')
                val: command parameter value (e.g. 'bbmapskimmer.sh')
        """
        pass

    def align(self, skip):
        """Performs the actual sequence alignment.

        Arguments:
            skip: Boolean - Whether or not to skip this particular alignment.
        """

        cline = shlex.split(self.command_line)
        wrapper = textwrap.TextWrapper(break_long_words=False,
                initial_indent = ' '*8,
                subsequent_indent = ' '*8)

        if VERBOSE:
            sys.stdout.write("\n")
            logging.info("Alignments...")
            if skip:
                pmsg("    Skipping the alignment. Would have called the following:")
            else:
                pmsg("    Running the alignment:")
            pmsg(wrapper.fill("{}".format(self.command_line)))

        if not skip:
            if (subprocess.call(cline) > 0):
                raise ExecutionException, \
                        ("Execution of the following command failed:\n{}".
                                format(self.command_line))

    def hits(self):
        """Generator that walks through each alignment hit.

        Implemented as a generator so as to avoid loading the entire alignment
        file into memory."""

        num_lines = sum(1 for line in open(self.results_file, 'rU') if line[0] != '@')
        indent = 8
        progress = Progress("Reading BBMap hits from file",
            10000, num_lines, indent=indent)
        with open(self.results_file,'rU') as f:
            for line in f:
                if line[0] == '@':
                    continue
                progress.progress()
                try:
                    # e.g. if self's class is BBMapExecution, then hit_cls will
                    # be BBmapExecutionHit
                    hit_cls_name = self.__class__.__name__ + 'Hit'
                    hit_cls = str_to_cls('prost.alignment', hit_cls_name)
                    hit = hit_cls(line)
                except UnmappedAlignmentException as e:
                    # Found an unmapped read.  Continue without creating the
                    # hit.
                    # TODO: If this is BBMap or BBMapSkimmer, this is likely a
                    # BBMap/BBMapSkimmer bug.  Report this bug.
                    perr(e)
                    continue

                # Filter out hits that have too many a) 3p mismatches,
                # b) non_3p mismatches, c) too many indel nts, d) any indel nts
                # if indels are not allowed.
                num_non_3p_mismatches = (hit.num_non_3p_mismatches
                        + (hit.num_indelnts * self.indelnt_penalty_multiplier))

                if (    hit.is_soft_clipped or
                        num_non_3p_mismatches > self.max_non_3p_mismatches or
                        hit.num_3p_mismatches > self.max_3p_mismatches or
                        (hit.has_indelnts and not self._allow_indels)):
                    hit.is_no_hit = True
                else:
                    hit.is_no_hit = False
                yield hit
        progress.done()

    def single_sequence_hit_sets(self):
        """Generator that walks through groups of search hits comprised of the
        same single query sequence.

        For example, assuming the following six search results:
            AAA ..., AAA ..., AAA ..., TGT ..., CCC ..., CCC ...
        Then this generator will yield three results:
            [AAA ..., AAA ..., AAA ...]
            [TGT ...]
            [CCC ..., CCC ...]

        Yields:
            [AlignmentExecutionHit]: A list of alignment hits which all contain
                the same query sequence.
        """

        hit_set = []
        last_query_sequence = None

        for hit in self.hits():
            if hit.query_sequence == last_query_sequence:
                hit_set.append(hit)
            elif last_query_sequence == None:
                hit_set.append(hit)
            else:
                # Current hit's query sequence is different from those of the
                # current set. Yield then begin building the next set.
                yield hit_set
                hit_set[:] = [hit]
            last_query_sequence = hit.query_sequence

        # Yield the last hit_set
        yield hit_set


class BBMapAlignmentExecution(AlignmentExecution):
    """A BBMap search.

    Example:
        bbmapskimmer.sh mdtag=t in="search.fa" maxindel=2
        out="alignment_results_bbmap_genome" slow=t outputunmapped=f idtag=f
        minid=0.50 ssao=f strictmaxindel=t usemodulo=f cigar=t sssr=0.25
        threads=12 path="/some/path/to/a/BBMap/DB" trimreaddescriptions=t
        secondary=t ambiguous=all maxsites=4000000 k=7 usejni=f
        maxsites2=4000000 idfilter=0.50
    """

    def generate_command_line(self):
        """Concrete implementation of AlignmentExecution.generate_command_line()."""

        # The command line string we will build up.
        self._cline = [EXECUTABLE_BBMAP]

        # First define the default parameters.
        params = {
            "mdtag": BBMAP_OPT_DEFAULT_MDTAG,
            "cigar": BBMAP_OPT_DEFAULT_CIGAR,
            "scoretag": BBMAP_OPT_DEFAULT_SCORETAG,
            "stoptag": BBMAP_OPT_DEFAULT_STOPTAG,
            "idtag": BBMAP_OPT_DEFAULT_IDTAG,
            "inserttag": BBMAP_OPT_DEFAULT_INSERTTAG,
            "usejni": BBMAP_OPT_DEFAULT_USEJNI,
            "outputunmapped": BBMAP_OPT_DEFAULT_OUTPUTUNMAPPED,
            "trimreaddescriptions": BBMAP_OPT_DEFAULT_TRIMREADDESCRIPTIONS,
            "k": BBMAP_OPT_DEFAULT_K,
            "usemodulo": BBMAP_OPT_DEFAULT_USEMODULO,
            "strictmaxindel": BBMAP_OPT_DEFAULT_STRICTMAXINDEL,
            "ambiguous": BBMAP_OPT_DEFAULT_AMBIGUOUS,
            "secondary": BBMAP_OPT_DEFAULT_SECONDARY,
            "sssr": BBMAP_OPT_DEFAULT_SSSR,
            "ssao": BBMAP_OPT_DEFAULT_SSAO,
            "maxsites": BBMAP_OPT_DEFAULT_MAXSITES,
            "maxsites2": BBMAP_OPT_DEFAULT_MAXSITES2,
            "slow": BBMAP_OPT_DEFAULT_SLOW,
        }

        # Set configurable params
        max_mismatches = max(self.max_3p_mismatches,
                self.max_non_3p_mismatches)
        params['out'] = self.results_file
        params['in'] = self.input_file
        params['threads'] = str(self._max_threads)
        params['path'] = os.path.expanduser(self.seq_db)
        if self._allow_indels:
            params['maxindel'] = str(
                self.max_non_3p_mismatches //
                self.indelnt_penalty_multiplier)
        else:
            params['maxindel'] = "0"


        # Calculate BBMap param 'minid' ("approximate minimum alignment
        # identity") and idfilter (independant of minid; sets exact minimum
        # identity allowed for alignments to be printed.) based upon the
        # configured max_mismatches and min_seq_length.  Hits longer than that
        # will be filtered out by prost.
        #
        # Actually, hard code this puppy for now.
        #params['minid'] = str((
        #    (  self._min_seq_length - self._max_mismatches) /
        #        float(self._min_seq_length)))
        params['minid'] = TODO_HARDCODED_MINID
        params['idfilter'] = params['minid']

        # Fourth, generate the command line string.
        for key in params.iterkeys():
            if(key in ("path", "in", "out")):
                # Put quotes around the paths of path, in, and out, which may
                # potentially have spaces in them
                self._cline.append('{}="{}"'.format(key, params[key].strip('"').strip()))
            else:
                self._cline.append('{}={}'.format(key, params[key].strip('"').strip()))
        self._cline = ' '.join(self._cline)


class AlignmentExecutionHit(SlotPickleMixin):
    """
    Abstract Base Class used to represent two concepts: a) an alignment (i.e.
    "hit") of a query sequence to a reference sequence, and b) the genomic
    location(s) of a ShortSeq.

    When used in the context of an alignment, the AlignmentExecutionHit is used
    to store, query, and manipulate information about the alignment (for
    example, see has_100_percent_core_identity()).

    When used in the context of a ShortSeq's genomic location(s), the
    name and location (i.e. genomic coordinates) of the AlignmentExecutionHit's
    reference sequence is defined as one of the ShortSeq's genomic locations.

    The only constraint is you need to use __slots__ for all
    ancestor/descendant classes.

    Attributes:
        designation: The Prost designation of this hit.
        is_no_hit: True if Prost post-filtered this hit (for example, if there
            were too many mismatches on the 3-prime end), False otherwise.

    """
    __metaclass__ = ABCMeta
    __slots__ = ('designation', 'is_no_hit')

    def __lt__(self, oth):
        """Less Than comparison function. Called by comparison operators (e.g. for
        sorting).  Compares the reference sequence.  Used to sort a ShortSeq's
        genomic_locations.

        Note:
            We want '+' strand hits to come before '-' strands hit (this is
            just an arbitrarily chosen convention).  To do so, we sort with
            "on_minus_strand", which places '+' strand first since False < True
            in python.

        Arguments:
            oth (AlignmentExecutionHit): The hit against which we are comparing
                this (self) hit.
        """
        return ((   self.reference_sequence_name,
                    self.on_minus_strand,
                    self.reference_start,
                    self.reference_end)
            <
                (   oth.reference_sequence_name,
                    oth.on_minus_strand,
                    oth.reference_start,
                    oth.reference_end))


    def __len__(self):
        """The 'length' of a short sequence's Genomic Location."""
        return 1 + abs(self.reference_start - self.reference_end)

    def __init__(self):
        self.designation = None

        # Ensure that additional common properties are set correctly.
        assert(type(self.reference_sequence_name) == str)
        assert(type(self.reference_start) == int)
        assert(type(self.reference_end) == int)
        assert(type(self.on_minus_strand) == bool)

    @abstractproperty
    def query_sequence(self):
        """Getter for the full original query sequence."""
        pass

    @abstractproperty
    def has_100_percent_core_identity(self):
        """Does this hit have 100% identity of its core sequence?

        Example:
            alignment   hit     ident
            AACCGGTTT   CCGG    100%
            TTCCGGCCC

        Only a mismatch on either/both ends and no mismatches in the middle.
        """
        pass

    @abstractproperty
    def is_full_length(self):
        """Is this hit a full length hit?

        Example:
            alignment   hit         ident    full_length?
            AACCGGTTT   AACCGGTTT   ! 100%   yes
            AACCTGTTT

        No mismatches on the ends.  Mismatches allowed in the middle. No
        indels.
        """
        pass

    @abstractproperty
    def has_indelnts(self):
        """Does this hit have any inserted or deleted nucleotides?"""
        pass

    @abstractproperty
    def num_indelnts(self):
        """Returns the number of inserted or deleted nucleotides in this hit.

        Example:
            20=2D8=3I20=  -> then num_indelnts = 2+3 = 5
        """
        pass

    @abstractproperty
    def num_3p_mismatches(self):
        """Returns the number of mismatches on the 3-prime side in this hit."""
        pass

    @abstractproperty
    def num_non_3p_mismatches(self):
        """Returns the number of mismatches that are not on the 3-prime side
        (i.e. those that are within the core sequence as well as on the 5-prime
        side).

        Note that this number does NOT consider indels.
        """
        pass

    @abstractproperty
    def alignment_length(self):
        """Returns the length of the alignment in this hit."""
        pass

    @abstractproperty
    def cigar5pto3p(self):
        """Return the CIGAR string reading from the 5p to 3p direction (i.e.
        strand-independent)."""
        pass

    @abstractproperty
    def has_5p_mismatch(self):
        """Does the query (i.e. short) sequence have a mismatch on its 5' end
        from the reference sequence?"""
        pass

    @abstractproperty
    def has_3p_mismatch(self):
        """Does the query (i.e. short) sequence have a mismatch on its 3' end
        from the reference sequence?"""
        pass

    @abstractproperty
    def reference_sequence_name(self):
        """Returns the reference sequence name (e.g. a linkage group) of this
        hit."""
        pass

    @abstractproperty
    def reference_start(self):
        """The 5-prime start on the reference sequence for this hit.  If the
        alignment is to the minus-strand of the reference, ref_start will be
        greater than ref_end."""
        pass

    @abstractproperty
    def reference_end(self):
        """The 3-prime start on the reference sequence for this hit.  If the
        alignment is to the minus-strand of the reference, ref_start will be
        greater than ref_end."""
        pass

    @abstractproperty
    def reference_start_with_clips(self):
        """The 5-prime start on the reference sequence for this hit adjusted
        for any 5-prime cliped nucleotides.  Used by (for example)
        is_seed_shifted(). See also reference_start above.  If the alignment is
        to the minus-strand of the reference, ref_start will be greater than
        ref_end.

        Examples:

            Example1::

                123456789012345678
                NNNCCAAATTGGCNNNNN   plus-strand
                    aaAAATTGGCaaa

                reference_start_with_clips =  5
                reference_start            =  7

            Example2::

                123456789012345678
                NNNCCAAATTGGCNNNNN   minus-strand
                    aaAAATTGGCaaa

                reference_start_with_clips =  17
                reference_start            =  14
        """
        pass

    @abstractproperty
    def reference_end_with_clips(self):
        """The 3-prime end on the reference sequence for this hit adjusted
        for any 3-prime cliped nucleotides.  Used by (for example)
        is_seed_shifted(). See also reference_start above.  If the alignment is
        to the minus-strand of the reference, ref_start will be greater than
        ref_end.

        Examples:

            Example::

                12345678901234567
                NNNCCAAATTGGCANNN   minus-strand
                aaAAATTGGCA

                reference_end_with_clips = 4
                reference_end            = 6
        """
        pass

    ## Genomic Location Properties ##

    @abstractproperty
    def lg(self):
        """Courtesy alias for reference_sequence_name().  Used by code to track
        the genomic locations of our short sequences."""
        pass

    @abstractproperty
    def start(self):
        """Courtesy alias for reference_sequence_start().  Used by code to track
        the genomic locations of our short sequences."""
        pass

    @abstractproperty
    def end(self):
        """Courtesy alias for reference_sequence_end().  Used by code to track
        the genomic locations of our short sequences."""
        pass

    ## Methods ##

    @abstractmethod
    def remove_query_sequence(self):
        """Remove the query_sequence to save memory."""
        pass

    def calculate_designation_full(self, designation_integer):
        """Calculates the full designation for sequences that have a fractional
        designation."""
        designation_fraction = \
            self.calculate_designation_fraction(designation_integer)
        return str(designation_integer) + '.' + str(designation_fraction)

    def calculate_designation_fraction(self, designation_integer):
        """For designations 3 and above, determine the fractional part of the
        sequence designation.

        For example, a designation '4' may become a designation '4.2'.

        Arguments:
            designation_integer: The part of the sequence designation to the
                left of the decimal point (e.g. for designation '3.2', the
                designation_integer is '3').
        Returns:
            str: The designation fraction, which is the part of the sequence
                designation to the right of the decimal point.
        """

        designation_fraction = ''
        if(self.has_5p_mismatch and self.has_3p_mismatch):
            designation_fraction = "5p3p"
        elif(self.has_5p_mismatch):
            designation_fraction = "5p"
        elif(self.has_3p_mismatch):
            designation_fraction = "3p"
        elif(self.is_full_length):
            designation_fraction = "0"
        if(self.has_indelnts):
            if designation_fraction == '':
                designation_fraction = INDEL_INDICATOR
            else:
                designation_fraction += ".{}".format(INDEL_INDICATOR)
        assert designation_fraction != None, \
            "unable to compute fractional part of sequence designation"
        return designation_fraction

    ## Genomic Location Methods ##

    def within_mirror(self, other, overlap=DEFAULT_MIRROR_OVERLAP):
        """Is this genomic location a potential mirror-mir of the 'other'
        genomic location?
        """
        within = False

        # same strand?
        if (    (self.start < self.end and other.start < other.end) or
                (self.start > self.end and other.start > other.end)):
            return False

        s_coords = sorted([self.start, self.end])
        o_coords = sorted([other.start, other.end])

        # Shrink the coords of 'other' by 'overlap' - 1.
        o_coords[0] += (overlap - 1)
        o_coords[1] -= (overlap - 1)

        # Aliases
        s_start = s_coords[0]
        o_start = o_coords[0]
        s_end = s_coords[1]
        o_end = o_coords[1]

        if self.lg != other.lg:
            # Have to be on the same lg to be a mirror.
            within = False
        elif overlap >= len(self) or overlap >= len(other):
            # If we ask for an overlap of 30 nts, but one of the mirs is only
            # 20nt, then we obviously can't overlap.
            within = False
        elif s_start >= o_start and s_start <= o_end:
            within = True
        elif s_end >= o_start and s_end <= o_end:
            within = True
        elif o_start >= s_start and o_start <= s_end:
            within = True
        elif o_end >= s_start and o_end <= s_end:
            within = True

        return within

    def within_wiggle(self, other, wiggle):
        """Is this genomic location within genetic wiggle of the 'other'
        genomic location?

        See 'ShortSeq.within_wiggle()' for more information.

        WARNING: This method is NOT symmetric.  If 'self' is within_wiggle of
        'other', it is not necessarily true that 'other' is within_wiggle of
        'self.  If 'self' is contained entirely within 'other', AND its ends
        are not within the wiggle length of 'other's ends, then self is
        within_wiggle of 'other', but the converse is not true.

        Arguments:
            other: GenomicLocation - A genomic location.
            wiggle: Integer - The wiggle room in base pairs.
        """

        within = False

        if self.lg != other.lg:
            within = False
        elif (abs(other.start - self.start) <= wiggle
                and abs(other.end - self.end) <= wiggle):
            # regular wiggle test
            within = True
        elif(other.start < other.end
                and self.start < self.end
                and other.start <= self.start
                and other.end >= self.end):
            # other 5'->3', self is contained entirely within other
            within = True
        elif(other.start > other.end
                and self.start > self.end
                and other.start >= self.start
                and other.end <= self.end):
            # other 3'->5', self is contained entirely within other
            within = True

        return within


class SamAlignmentExecutionHit(AlignmentExecutionHit):
    """Concrete implementation of AlignmentExecutionHit, SAM file version."""

    #
    #               === Handy Notes about SAM files ===
    #
    #
    # Note: Copied from unlicensed version of SAMv1.tex found at: https://github.com/samtools/hts-specs/blob/90b9d9edad9cd9c7c07318dbb37cacc543d0e62e/SAMv1.tex
    #
    #
    # The example from the SAM Format Specification (v1.4?, 28 Dec 2014):
    #
    # Suppose we have the following alignment with bases in lower cases clipped from
    # the alignment. Read r001/1 and r001/2 constitute a read pair; r003 is a
    # chimeric read; r004 represents a split alignment.
    #
    # Coor      12345678901234 5678901234567890123456789012345
    # ref       AGCATGTTAGATAA**GATAGCTGTGCTAGTAGGCAGTCAGCGCCAT
    # +r001/1         TTAGATAAAGGATA*CTG
    # +r002          aaaAGATAA*GGATA
    # +r003        gcctaAGCTAA
    # +r004                      ATAGCT..............TCAGC
    # -r003                             ttagctTAGGC
    # -r001/2                                         CAGCGGCAT
    #
    # The corresponding SAM format is:
    #
    # @HD VN:1.5 SO:coordinate
    # @SQ SN:ref LN:45
    # r001   99 ref  7 30 8M2I4M1D3M = 37  39 TTAGATAAAGGATACTG *
    # r002    0 ref  9 30 3S6M1P1I4M *  0   0 AAAAGATAAGGATA    *
    # r003    0 ref  9 30 5S6M       *  0   0 GCCTAAGCTAA       * SA:Z:ref,29,-,6H5M,17,0;
    # r004    0 ref 16 30 6M14N5M    *  0   0 ATAGCTTCAGC       *
    # r003 2064 ref 29 17 6H5M       *  0   0 TAGGC             * SA:Z:ref,9,+,5S6M,30,1;
    # r001  147 ref 37 30 9M         =  7 -39 CAGCGGCAT         * NM:i:1

    # Sam Mandatory Fields
    # Col Field Type   Regexp/Range            Brief description
    # 1   QNAME String [!-?A-~]{1,255}         Query template NAME
    # 2   FLAG  Int    [0,2^16-1]              bitwise FLAG
    # 3   RNAME String \*|[!-()+-<>-~][!-~]*   Reference sequence NAME
    # 4   POS   Int    [0,2^31-1]              1-based leftmost mapping POSition
    # 5   MAPQ  Int    [0,2^8-1]               MAPping Quality
    # 6   CIGAR String \*|([0-9]+[MIDNSHPX=])+ CIGAR string
    # 7   RNEXT String \*|=|[!-()+-<>-~][!-~]* Ref. name of the mate/next read
    # 8   PNEXT Int    [0,231-1]               Position of the mate/next read
    # 9   TLEN  Int    [-2^31+1,2^31-1]        observed Template LENgth
    # 10  SEQ   String \*|[A-Za-z=.]+ segment  SEQuence
    # 11  QUAL  String [!-~]+                  ASCII of Phred-scaled base QUALity+33

    # FLAG field  Bit Description
    # 0x1         template having multiple segments in sequencing
    # 0x2         each segment properly aligned according to the aligner
    # 0x4         segment unmapped
    # 0x8         next segment in the template unmapped
    # 0x10        SEQ being reverse complemented
    # 0x20        SEQ of the next segment in the template being reversed
    # 0x40        the first segment in the template
    # 0x80        the last segment in the template
    # 0x100       secondary alignment
    # 0x200       not passing quality controls
    # 0x400       PCR or optical duplicate
    # 0x800       supplementary alignment

    # CIGAR String Operations
    # Op BAM Description
    # M   0  alignment match (can be a sequence match or mismatch)
    # I   1  insertion to the reference
    # D   2  deletion from the reference
    # N   3  skipped region from the reference
    # S   4  soft clipping (clipped sequences present in SEQ)
    # H   5  hard clipping (clipped sequences NOT present in SEQ)
    # P   6  padding (silent deletion from padded reference)
    # =   7  sequence match
    # X   8  sequence mismatch
    #
    # Notes on CIGAR String Operations
    # * H can only be present as the first and/or last operation.
    # * S may only have H operations between them and the ends of the CIGAR string.
    # * For mRNA-to-genome alignment, an N operation represents an intron. For
    #   other types of alignments, the interpretation of N is not defined.
    # * Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.

    # MD TAG Specification
    # Tag Type Description
    # MD    Z  String for mismatching positions. Regex : [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*6
    #
    # Notes on MD TAG Specification
    # The MD field aims to achieve SNP/indel calling without looking at the
    # reference. For example, a string '10A5^AC6' means from the leftmost
    # reference base in the alignment, there are 10 matches followed by an A on
    # the reference which is different from the aligned read base; the next 5
    # reference bases are matches followed by a 2bp deletion from the
    # reference; the deleted sequence is AC; the last 6 bases are matches. The
    # MD field ought to match the CIGAR string.


    # Defines the required and optional SAM file fields.
    __sam_fields = tuple('qname flag rname pos mapq cigar rnext pnext tlen seq qual tags'.split())
    __sam_fields_retained = ('qname', 'flag', 'rname', 'pos', 'cigar')
    __additional_properties = ('mdtag', 'cigar_tokens', 'is_soft_clipped')
    # __slots__ = __sam_fields + __additional_properties # uncomment for debugging
    __slots__ = __sam_fields_retained + __additional_properties

    ## API properties ##

    @property
    def query_sequence(self):
        return self.qname

    @property
    def reference_sequence_name(self):
        return self.rname
    # Alias
    lg = reference_sequence_name

    # Example of determining the rightmost reference position in an alignment:
    # left_pos  = Given by SAM alignment (POS field)
    # right_pos = left_pos + ((sum of numbers before [M=XDN]) - 1
    #             Note that we ignore [ISHP] as do not affect the current
    #             position in the reference.
    #
    # 123456789--0123456789012345
    # NNNCCAAAT--TGGCATGCANNNNN                 2S4M2I5M2D2M  2S4=2I2=1X2=2D2=
    #    aaAAATCGTGACA--CA                      0C5A2^TG1
    #
    # MD Tag:    0C5A2^TG1
    # old cigar: 2S4M2I5M2D2M
    # ext cigar: 2S4=2I2=1X2=2D2=
    #
    # left_pos = 6
    # right_pos (old cigar) = 6 + (4 + 5 + 2 + 2) - 1 = 18
    # right_pos (ext cigar) = 6 + (1 + 3 + 2 + 1 + 2 + 2 + 2) - 1 = 18

    @property
    def reference_start(self):
        # plus-strand: POS
        # minus-strand: POS + (sum of numbers before [M=XDN]) - 1
        # See example above.
        start = None
        if self.on_minus_strand:
            tokens = self.cigar_tokens
            nts = (t[0] for t in tokens if t[1] in ('M', '=', 'X', 'D', 'N'))
            read_len = sum(nts)
            start = self.pos + read_len - 1
        else:
            start = self.pos

        return start
    # Alias
    start = reference_start

    @property
    def reference_end(self):
        # plus-strand: POS + (sum of numbers before [M=XDN]) - 1
        # minus-strand: POS
        # See example above.
        start = None
        if not self.on_minus_strand:
            tokens = self.cigar_tokens
            nts = (t[0] for t in tokens if t[1] in ('M', '=', 'X', 'D', 'N'))
            read_len = sum(nts)
            start = self.pos + read_len - 1
        else:
            start = self.pos

        return start
    # Alias
    end = reference_end

    @property
    def reference_start_with_clips(self):
        # Used only for new categories (e.g. %_seed_edited)
        # plus-strand: POS - (sum of 5prime [SH]s)
        # minus-strand: POS + (sum of numbers before [M=XDN]) + (sum of 5prime [SH]s) - 1
        start = None
        tokens = self.cigar_tokens
        if self.on_minus_strand:
            # Go in at most two tokens to find only 5' clipped nts
            clipped_5p_H = tokens[-1][0] if tokens[-1][1] == 'H' else 0
            idx_S = -1 if clipped_5p_H == 0 else -2
            clipped_5p_S = tokens[idx_S][0] if tokens[idx_S][1] == 'S' else 0
            clipped_5p_nts = clipped_5p_H + clipped_5p_S
        else:
            # Go in at most two tokens to find only 5' clipped nts
            clipped_5p_H = tokens[0][0] if tokens[0][1] == 'H' else 0
            idx_S = 0 if clipped_5p_H == 0 else 1
            clipped_5p_S = tokens[idx_S][0] if tokens[idx_S][1] == 'S' else 0
            clipped_5p_nts = clipped_5p_H + clipped_5p_S

        if self.on_minus_strand:
            nts = (t[0] for t in tokens if t[1] in ('M', '=', 'X', 'D', 'N'))
            align_len = sum(nts)
            start = self.pos + align_len - 1 + clipped_5p_nts
        else:
            start = self.pos - clipped_5p_nts

        return start

    @property
    def reference_end_with_clips(self):
        query_seq_len = sum([t[0] for t in self.cigar_tokens])

        if self.on_minus_strand:
            ref_end = self.reference_start_with_clips - query_seq_len + 1
        else:
            ref_end = self.reference_start_with_clips + query_seq_len - 1

        if ref_end <= 0:
            # Hopefully should never happen
            raise Exception

        return ref_end

    @property
    def has_indelnts(self):
        # True if there is an I or D present.
        i_and_ds = [t[1] for t in self.cigar_tokens if t[1] in ('I', 'D')]
        return True if len(i_and_ds) else False

    @property
    def num_indelnts(self):
        # Count the number of 'I's and 'D's. - Note that mismatches are *not*
        # counted.
        i_and_ds = [t[0] for t in self.cigar_tokens if t[1] in ('I', 'D')]
        return (sum(i_and_ds))

    @property
    def alignment_length(self):
        # sum of numbers before [M=XDN]
        tokens = self.cigar_tokens
        return sum(t[0] for t in tokens if t[1] in ('M', '=', 'X', 'D', 'N'))

    @property
    def cigar_tokens_5pto3p(self):
        tokens = self.cigar_tokens
        if self.on_minus_strand:
            tokens = tokens[::-1]
        return tokens

    @property
    def cigar5pto3p(self):
        tokens = self.cigar_tokens
        if self.on_minus_strand:
            tokens = tokens[::-1]
        # flatten and join into one big string
        return "".join(str(i) for i in itertools.chain.from_iterable(tokens))

    ## Internal properties ##

    @property
    def _flag_unmapped(self):
        """Does this hit have the SAM flag "0x4 segment unmapped" set?"""
        return 0x004 == (0x004 & self.flag)

    @property
    def _flag_reverse_complement(self):
        """Does this hit have the SAM flag "0x10 SEQ being reverse
        complemented" set?"""
        return 0x010 == (0x010 & self.flag)

    @property
    def _flag_secondary_alignment(self):
        """Does this hit have the SAM flag "0x100 secondary alignment" set?"""
        return 0x100 == (0x100 & self.flag)

    @property
    def on_minus_strand(self):
        """Is this hit on the minus strand?"""
        return 0x010 == (0x010 & self.flag)
    def __init__(self, hit_line):

        # MDTags are not always present.  Set to ""
        self.mdtag = ''

        # set self's attributes to those found in the SAM hit
        aligned_num_fields = 12 # should be a constant...
        vals = hit_line.rstrip('\n').split('\t', aligned_num_fields - 1)
        if len(vals) == (aligned_num_fields - 1):
            # When it doesn't align, then tags field is empty.
            (self.qname, self.flag, self.rname, self.pos, self.cigar) = \
                    operator.itemgetter(0, 1, 2, 3, 5)(vals)
            self.mdtag = ""
        else:
            (self.qname, self.flag, self.rname, self.pos, self.cigar, tags) = \
                    operator.itemgetter(0, 1, 2, 3, 5, 11)(vals)
            for tag in tags.split('\t'):
                if "MD:Z:" == tag[0:5]:
                    self.mdtag = tag[5:]
                    break
        # Make ints ints
        self.flag = int(self.flag, 0)
        self.pos = int(self.pos, 0)

        # # Uncomment to collect everything else for debugging purposes
        # (self.mapq, self.rnext, self.pnext, self.tlen, self.seq, self.qual) = \
        #        operator.itemgetter(4, 6, 7, 8, 9, 10)(vals)

        # Ensure that we have an alignment. See hits() where this exception is
        # caught.
        #
        # Originally added because BBMapSkimmer was outputting a single flawed
        # alignment:
        #   TACTCAGGATCTAGTTTTCCAACTTTGC 272 17 50130514 24 * * 0 0 * * AM:i:24 NH:i:2 YR:i:2090
        if self.cigar == '*':
            raise UnmappedAlignmentException, \
                """Warning: Found an unmapped alignment in your SAM file.
                    Either this is a bug in your aligner, or you should
                    configure your aligner to not output unmapped reads. The
                    offending line is '{}'.""".format(hit_line.rstrip('\n'))

        # Unlikely, but blow up if you find any of these flags
        for flag in (0x1, 0x2, 0x8, 0x20, 0x40, 0x80, 0x200, 0x400, 0x800):
            if flag == (flag & self.flag):
                raise Exception, \
                    "ERR925: Unknown SAM Flag {} in hit:\n{}".format(
                        flag, hit_line)

        # Memoize cigar_tokens (about 3.2 times faster than w/o memoizing)
        self.cigar_tokens = self.parse_cigar_tokens()

        ## Set common properties
        self.is_soft_clipped = any(t[1] == 'S' for t in self.cigar_tokens)

        # Run super to ensure certain requirements are met (e.g.
        # 'on_minus_strand' is set)
        super(SamAlignmentExecutionHit, self).__init__()

    def __repr__(self):
        strand = "-" if self.on_minus_strand else "+"
        return '{}:{}-{} ({}) | {} {} ({})'.format(
                self.lg, self.start, self.end, strand, self.cigar5pto3p,
                self.mdtag, self.query_sequence)

    def remove_query_sequence(self):
        self.qname = None

    def parse_cigar_tokens(self):
        """Parse the cigar string into tokens.

        Returns:
            ((int,String),): A tuple of (int, String) tuples, where each
                (int, String) tuple represents one cigar token. For example, a
                cigar string of "17=1X5=" will be parsed into
                ((17, '='), (1, 'X'), (5, '=')).
        """
        # cigar_str = "10=1X9="
        # l = ['10', '=', '1', 'X', '9', '=', '']
        l = re.split('([MIDNSHP=X])', self.cigar)
        # Yup, this works.
        return tuple((int(n),s) for n,s in itertools.izip(l[::2], l[1::2]))

    def assumption_assertions_cigar_strings(self):
        """
        We make several assumptions about the characteristics of
        cigar strings.  These assumptions are enforced here.
        """
        report = "\tQuery:\t{}\n\tCIGAR:\t{}\n".format(self.qname, self.cigar)

        # If it doesn't hit, it doesn't have a cigar string.
        if self.cigar == '*':
            return

        # Assert no 'I' or 'D' at beginning or end of ext cigar str.
        # This should presumably be impossible.
        for idx in (0, -1):
            if self.cigar_tokens[idx] in ('I', 'D'):
                raise SamExtendedCigarParsingException, \
                    """ERR503: Oops, you found an unexpected 'I' or 'D' at the
                        beginning or end of the extended cigar string, which
                        Prost does not currently support, as this doesn't
                        really make sense. Please contact the Prost developers
                        with this report (copy and paste it) so that we may
                        examine it. Thank you.\n{}""".format(report)

        # Assert no 'N' or 'P' or 'H'.
        for t in self.cigar_tokens:
            if t[1] in (
                'N',   # We don't expect introns in this dataset.
                'P',   # TODO: Really? We don't expect padding in our dataset?
                'H'    # TODO: Seems reasonable, no chimeric reads in our dataset?
            ):
                raise SamExtendedCigarParsingException, \
                    """ERR503: Oops, you found a CIGAR operation ('N', 'P', 'H')
                        that Prost doesn't support. Please contact the
                        Prost developers and we'll add support for
                        this.\n{}""".format(report)

        # TODO: comment out this block after confirming that BBMap
        # parsed ext cigar strings do indeed match sequence length.

        # # Uncomment this assertion to test if an aligner's ext cigar str
        # # representation of sequence length matches the actual sequence length
        # # (as specified by length(self.qname)).
        # nts = (t[0] for t in self.cigar_tokens if t[1] in ('M', 'I', 'S', '=', 'X'))
        # read_len = sum(nts)
        # if read_len != len(self.qname):
        #   raise SamExtendedCigarParsingException, \
        #       "ERR503: TEMPORARY: len of parsed(cigar) != len(seq).\n{}".format(report)


class SamExtendedCigarAlignmentExecutionHit(SamAlignmentExecutionHit):
    """A SamAlignmentExecutionHit which internally parses a SAM extended CIGAR
    string."""

    __slots__ = ()

    @property
    def has_100_percent_core_identity(self):
        # True if there's only one '=' in the entire sequence.
        equals = [t[1] for t in self.cigar_tokens if t[1] == '=']
        return len(equals) == 1

    @property
    def is_full_length(self):
        # If the first or last cigar token is not an '=', it's not full length.
        # Do not allow indels.
        if self.cigar_tokens[0][1] != '=' or self.cigar_tokens[-1][1] != '=':
            return False
        elif self.has_indelnts:
            return False
        else:
            return True

    @property
    def _num_mismatches(self):
        # Count the number of MXSH's. - Note that indels are *not* counted.
        # Note: indels are not counted.
        # This is NOT a public method, don't use it.
        mxsh = [t[0] for t in self.cigar_tokens if t[1] in ('M', 'X', 'S', 'H')]
        return sum(mxsh)

    @property
    def num_3p_mismatches(self):
        """
        Note:
            We define mismatches and indels to be mutually exclusive.

            In the extended cigar format, we expect matches to be '=' and
            mismatches to be 'X'.  Because of this, we expect no 'M's (or at
            least, an exceedingly rare number).

            In our dataset, BBMap actually produces less than a handful of
            'M' cigar operations.  These occur when the read aligns to a
            segment in the reference which is on the edge of many 'N's.  For
            the portion that aligns to [ACGT], those nts show up as '=', for
            the portion that overlaps the 'N's, those nts show up as 'M' (e.g.
            '21M7=').  This is a little odd, but simply count them as a
            mismatch for now.  It's rare anyway.

            The remainder (i.e. [XSH]) are obviously mismatches.

            P is similar enough to an I or D in this case, so ignore all three
            of those ([PID]).
        """
        total = 0
        if self.on_minus_strand:
            ctokens = self.cigar_tokens
        else:
            ctokens = self.cigar_tokens[::-1]
        # Walk through cigar tokens 3p-to-5p until you hit '='; collect and add
        # all nts before that.
        for num, kind in ctokens:
            if kind == '=':
                break
            total += num
        return total

    @property
    def num_non_3p_mismatches(self):
        return self._num_mismatches - self.num_3p_mismatches

    @property
    def has_5p_mismatch(self):
        # plus-strand: False if it doesn't begin with '='.
        # minus-strand: False if it doesn't end with '='.
        # Note: __init__ assertion prevents leading/trailing I's and D's.
        idx = -1 if self.on_minus_strand else +0
        return self.cigar_tokens[idx][1] != '='

    @property
    def has_3p_mismatch(self):
        # plus-strand: False if it doesn't end  with '='.
        # minus-strand: False if it doesn't begin with '='.
        # Note: __init__ assertion prevents leading/trailing I's and D's.
        idx = +0 if self.on_minus_strand else -1
        return self.cigar_tokens[idx][1] != '='


class BBMapAlignmentExecutionHit(SamExtendedCigarAlignmentExecutionHit):
    __slots__ = ()

    def __init__(self, hit_line):
        super(BBMapAlignmentExecutionHit, self).__init__(hit_line)
        self.assumption_assertions_cigar_strings()


# Note - consider making SamMDTagHit inherit from SamExtCigHit, and then simply
#        convert from OldCig/MDTag to ExtCig (in init()), and then all other
#        code (e.g. is_full_length, has_3p_overhang or whatever) will be based
#        on parsing the ExtCig instead of trying to mirror such code by parsing
class SamMDTagAlignmentExecutionHit(SamExtendedCigarAlignmentExecutionHit):

    ## Internal properties ##

    @property
    def mdtag_tokens(self):
        """Parse the mdtag string into tokens.

        Returns:
            (int|String): A tuple of with elements of type int or String.  The
                As does the MD tag does in the SAM file, the mdtag_tokens tuple
                starts and ends with an int, and alternates between ints and
                Strings.  For example, an mdtag string of '14^TCCAC7A0A10' is
                parsed into (14, '^TCCAC', 7, 'A', 0, 'A', 10).

        """
        # mdtag_str = '14^TCCAC7AA10'
        # l = ['14', '^TCCAC', '7', 'AA', '10']
        l = re.split('([0-9]+)', mdtag)[1:-1]
        for i in xrange(0, len(l), 2):
            l[i] = int(l[i])
        # l = [14, '^TCCAC', 7, 'AA', 10]
        return tuple(l)

# vim: softtabstop=4:shiftwidth=4:expandtab
