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


"""
Prost!
    1. the usual toast when drinking alcohol; cheers
    2. PRocessor Of Small Transcripts

A script to quantify and annotate miRNA expression.
"""

# Python 3 imports
from __future__ import absolute_import
from __future__ import division

# version
from prost._version import __version__
__version_info__ = tuple(__version__.split('.'))

# Prost imports
from prost.constants import *
from prost.common import (
    ControlFlowException,
    ConfigurationException,
    CannotContinueException,
    PrerequisitesException,
    MirModificationCalculationException,
    ModificationThingEncounteredNAlignment,
    ArgumentTypeException,
    pmsg,
    perr,
    SlotPickleMixin,
    Progress)
from prost.alignment import (Alignments, GenomeAlignment)
from prost.annotation import (Annotation,
                              MirbaseAnnotation,
                              MirbaseMirAnnotation,
                              MirbaseMirReverseAnnotation,
                              MirbaseHairpinAnnotation,
                              BiomartOtherRNAAnnotation)
#from prost.db_caching import (
#    DB_PROXY,
#    SequenceDB,
#    ShortSeqCache,
#    GenomicLocationCache)
import prost.excel

# Other imports
import sys
import shlex
import subprocess
import collections

# For prost.log logging
import logging
import logging.config

# for debugging:
import pprint
# import json, yaml
import traceback

# for abstract base classes:
from abc import ABCMeta, abstractmethod, abstractproperty

# for command line arguments
import argparse

# For prost configuration file
import ConfigParser

# Regular expressions (for parsing config file)
import re

# For pickling
import cPickle as pickle
import os.path

# For file caching
from contextlib import contextmanager

# For more efficient SamplesCounts
from array import array

# For iterating!
import itertools

# For timing things
import time

# For quicker sorting (operator.[item|attr]getter)
import operator

# For log fold changes (e.g. in arm_switch_candidates)
import math

# For determining system memory size
import psutil

#################
### Functions ###
#################

def get_85perc_system_memory():
    """Returns an int representing roughly 85% of the system memory.

    Returns:
        int: 85% of the system memory.
    """
    return int(psutil.virtual_memory().total * 0.85 / 1024**3)

def setup_logging(conf, time_start):
    """Setup prost.log logging.

    Does the following:
        - configures the logging
        - adds a "header" line that includes the version and date
        - adds a line showing how Prost was executed from the command line
        - dumps contents of prost.config and samples_filelist to log

    Arguments:
        conf (Configuration): The Configuration singleton object.
        time_start (float): The time in seconds since the epoch.

    """
    # Configure logging
    logging.config.dictConfig(LOGGING)

    # First line of log, version + timedate.
    left = "Prost! version {}.".format(__version__)
    right_len = 79 - len(left)
    pmsg("{}{:>{}}.".format(
        left, time.strftime('%l:%M:%S %p %Z on %b %d, %Y'), right_len))
    pmsg("")

    # Next, add how Prost was run from the command line (only to the log).
    args = sys.argv
    args[0] = os.path.basename(sys.argv[0])
    logging.info("Excecuted as:\n")
    logging.info("    {}".format(" ".join(args)))

    # Next, add the full contents of the configuration file.
    header = "=== Configuration File '{}' Contents ===".format(
        conf.general.config_file)
    logging.info("")
    logging.info("{:^79}".format(header))
    logging.info("")
    with open(os.path.expanduser(conf.general.config_file), 'r') as f:
        logging.info(f.read())

    # Next, add the full contents of the samples_filelist.
    header = "=== Samples Filelist '{}' Contents ===".format(
        conf.general.samples_filelist)
    logging.info("")
    logging.info("{:^79}".format(header))
    logging.info("")
    with open(os.path.expanduser(conf.general.samples_filelist), 'r') as f:
        logging.info(f.read())

    # Finally, add a header line for Prost exection.
    header = "=== Prost Execution ==="
    logging.info("{:^79}".format(header))
    logging.info("")

def to_mirror_fingerprint(gen_loc, bucket_length=DEFAULT_MIRROR_BUCKET_LENGTH):
    """Creates a fingerprint to be used by mirror miR detection.

    Args:
        gen_loc (GenomicLocation): The list of genomic locations.
        bucket_length (int): The length in nucleotides of the buckets.

    Returns:
        ((str, int), (str, int)): A low and high mirror fingerprint for the
            genomic location passed. Each fingerprint is of type:

                (A, B)

            where types are:

                A: str
                B: int

            and where:

                'A' is the linkage groups of gen_locs,
                'B' is either "low" or "high" rounded start site,

    """

    # Do not allow bucket lengths to be odd.  Haven't thought much about this,
    # but want to make sure nothing odd happens.
    if ((bucket_length % 2) == 1):
        raise CannotContinueException, \
            """Cannot continue: MIRROR_BUCKET_LENGTH cannot be odd."""

    # Bucket length and bucket offset
    bl = bucket_length
    bo = bucket_length//2

    gl_low = (gen_loc.lg, gen_loc.start//bl*bl)
    gl_high = (gen_loc.lg, (gen_loc.start+bo)//bl*bl)
    return (gl_low, gl_high)

def to_wiggle_fingerprints(gen_locs, max_seq_len):
    """Builds two wiggle fingerprints (one "low", and one "high") based upon the
    genomic locations passed.  These fingerprints are used to bucket together
    short sequences that have a greater likelihood of being within_wiggle of one
    another.

    The two fingerprints are tuples of the following "subfingerprints":

    1. The starting basepair of the first genomic location.
    2. A sorted tuple of all the linkage groups of the genomic locations.
    3. A tuple of the strands on which the genomic locations are found (i.e.
        either 5'->3' or 3'->5').

    The 2nd and 3rd subfingerprint above are identical between the two
    fingerprints, while the 1st subfingerprint above differs.

    The linkage groups are split into bucket sizes of 4*max_seq_len.  The "low"
    fingerprint starts at 0 in a sequence like so (where 'm' is max_seq_len):

        0, 4m, 8m, 12m, ...

    For the first genomic location in gen_locs has 'start' < 4m, it is put in
    the first bucket; if 4m <= 'start' < 8m, it is put in the second bucket,
    and so on.

    The buckets corresponding to the "high" fingerprint are offset from those
    corresponding to the "low" buckets by 2*max_seq_len.  The sequence for the
    high fingerprint is like so:

        0, 2m, 6m, 10m, 14m, ...

    Example:
        Assume two genomic locations gl1 and gl2 and max_seq_length = 30, with
        gl1.start = 29 and gl2.start = 31.  Then their starting basepair "low"
        subfingerprints are 0 and 30 respectively, while their "high"
        fingerprints are 30 and 30 respectively.

    TODO: Could make this more memory efficient (and possibly the dumping to
        file cache more time efficient) by taking only the first three genomic
        locations to create the fingerprint (or even the first 10).  Really,
        once you get past several genomic locations, there's no point in storing
        all those in the key hash.

    Args:
        gen_locs ([GenomicLocation]): The list of genomic locations on which to
            calculate the two wiggle_fingerprints.
        max_seq_length (int): The maximum sequence length of a ShortSeq (as
            configured).

    Returns:
        ((int, (int,), (int,)), (int, (int,), (int,))): Two wiggle fingerprints,
            one "low", one "high".  A wiggle fingerprint is a tuple of type:

                (A, B, C)

            where types are:

                A: int
                B: tuple of ints
                C: tuple of ints

            and where:

                'A' is either "low" or "high" fingerprint,
                'B' is the linkage groups of gen_locs,
                'C' is the strandedness of the gen_locs.
    """

    assert (len(gen_locs) > 0), \
            """to_wiggle_fingerprints() requires a non-empty list of genomic
            locations."""


    # Bucket length and bucket offset
    bl = max_seq_len*4
    bo = max_seq_len*2

    # The only base pair location we fingerprint is the start of the first
    # genomic_location.
    gl0 = gen_locs[0]
    gl0_low =  gl0.start//bl*bl
    gl0_high = (gl0.start+bo)//bl*bl

    # Linkage group fingerprint
    lg_fingerprint = tuple(gl.lg for gl in gen_locs)

    # Strand fingerprint ('1' means 5'->3', '-1' means 3'->5')
    strand_fingerprint = \
            tuple(sorted(((1 if (gl.start - gl.end < 0) else -1) for gl in gen_locs)))

    # Build the two fingerprints
    fp_low = (gl0_low, lg_fingerprint, strand_fingerprint)
    fp_high = (gl0_high, lg_fingerprint, strand_fingerprint)
    return (fp_low, fp_high)

@contextmanager
def file_caching(conf, stage, *obj_list):
    """Context manager that provides file caching of objects at various stages.

    Currently pretty dumb, doesn't invalidate caches (if the files are present,
    it will use them).  Still, a nice time saver for the developer.  In the
    future, could rework it to date stamp the cache files.  Might be nice in
    the near future to have it *always* write cache files, and the command line
    flags would turn on/off whether those cache files are used.

    Arguments:
        conf: The Configuration singleton object.
        stage: The number of the stage we are caching.


    Usage is like so::

        with file_caching(_conf, 4, large_obj1, large_obj2) as cached_objs:
            if not cached_objs:
                large_obj1.expensive_task3()
                large_obj2.expensive_task4()
            else:
                large_obj1 = cached_objs[0]
                large_obj2 = cached_objs[1]

    The example above would produce two cache files:
        .file_cache_stage4_object0.p
        .file_cache_stage4_object1.p
    """

    file_name = ".file_cache_stage{}_object{}.p"
    cache_files_ok = True
    cached_objs = []

    # Check if pickle files are present
    for i in xrange(len(obj_list)):
        if not os.path.isfile(file_name.format(stage, i)):
            cache_files_ok = False

    if cache_files_ok and conf.general.read_from_file_cache:
        # Do not execute, instead, load from cache.
        progress = Progress("Caching: Loading stage {} from cache file(s)".
                format(stage))
        for i in xrange(len(obj_list)):
            cache_file = file_name.format(stage, i)
            cached_objs.append(pickle.load(open(cache_file, "rb")))
        yield cached_objs
        progress.done()
    elif conf.general.write_to_file_cache:
        # Regular execution, and dump to cache
        yield False
        progress = Progress("Caching: Dumping stage {} to cache file(s)".
                format(stage))
        for i, obj in enumerate(obj_list):
            cache_file = file_name.format(stage, i)
            pickle.dump(obj, open(cache_file, "wb"))
        progress.done()
    else:
        # Regular execution, no caching
        yield False


def annotations_of(annotation_cls, short_seqs, kind, *seq_str_list):
    """Returns a set() of annotations of the passed list of short sequences that
    that conform to the passed annotation class and annotation kind.

    Args:
        annotation_cls (type): The class of the annotation.
        short_seqs (ShortSeqs): The ShortSeqs singleton object.
        kind (Annotation.Kind): The kind of annotation (e.g.
            other_species_hairpin).  if kind is None, then the 'kind' filter is
            ignored.
        *seq_str_list ([str]): A list of ShortSeq sequence strings for which the
            set() of annotations is desired.

    Returns:
        set(): A set of annotations.

    Example:
        # Return the species hairpin annotations for the list of sequences.
        bn = some Bin
        short_seqs = The ShortSeqs singleton.
        anno_class = MirbaseHairpinAnnotation
        anno_kind = MirbaseHairpinAnnotation.Kind.species_hairpin
        annos = annotations_of(anno_class, short_seqs, anno_kind, \*(bn.joining_short_seq_strs))

    """

    annos = set()
    for short_seq in (short_seqs[s] for s in seq_str_list):
        for anno in short_seq.annotations:
            if (anno.__class__ == annotation_cls):
                if (kind):
                    if (anno.kind == kind):
                        annos.add(anno)
                else:
                    annos.add(anno)
    return annos

def log_fold_changes_expr_5p_3p(norms_5p, norms_3p):
    """Calculate per-sample log fold change of 5p/3p expressions.

    Notes:
        * If both expr_5p and expr_3p are 0, then log_fold_change = 0.
        * If expr_5p != 0 and expr_3p == 0, then log_fold_change = 1.
        * If expr_5p == 0 and expr_3p != 0, then log_fold_change = -1.

    Args:
        norms_5p (SamplesCountsWithNorms): The normalized per-sample counts of
            the 5p mature mir.
        norms_3p (SamplesCountsWithNorms): The normalized per-sample counts of
            the 3p mature mir.

    Returns:
        [floats]: An array of floats representing per-sample log fold change of
            5p/3p expression ranges.

    """
    if (type(norms_5p) != SamplesCountsWithNorms or
            type(norms_3p) != SamplesCountsWithNorms):
        raise ArgumentTypeException

    log_fold_changes = []
    for i in xrange(norms_5p.num_samples):
        expr_5p = norms_5p._norms[i]
        expr_3p = norms_3p._norms[i]
        if expr_5p == 0 and expr_3p == 0:
            change = 0
        else:
            if expr_5p == 0:
                change = -1
            elif expr_3p == 0:
                change = 1
            else:
                change = math.log(expr_5p/expr_3p, 2)
        log_fold_changes.append(change)

    return log_fold_changes


##############################################
### Config File and Command Line Arguments ###
##############################################

class Configuration(object):

    _defaults = { }

    @property
    def annotation_alignment_names(self):
        """Returns the section names of the annotation alignments."""
        annos = []
        for section_name in self._config.sections():
            if re.search(r'^AnnotationAlignment', section_name):
                annos.append(section_name)
        return annos

    @property
    def alignments(self):
        """Returns a dict representing configured options in the
        'GenomeAlignment' and 'AnnotationAlignment[0-9][0-9] sections of the
        config_file.

        The key in 'alignments' is the name of the alignment section (e.g.
        'AnnotationAlignment2'), and the value is another dict representing the
        configured options in that alignment section.
        """

        aligns = collections.OrderedDict()
        aligns['GenomeAlignment'] = self.get_section('GenomeAlignment')
        for section in self._config.sections():
            if re.search(r'^AnnotationAlignment', section):
                aligns[section] = self.get_section(section)
        return aligns

    def __init__(self):
        self._argparser = argparse.ArgumentParser(
                description=__doc__,
                formatter_class=argparse.RawDescriptionHelpFormatter,
            )
        # Note that we allow "no_value" fields here only for cleaner error
        # messages.  Later, in get_section(), we explicitly forbade "no_value"
        # fields.
        self._config = ConfigParser.SafeConfigParser(allow_no_value=True)

    def _setup_argparser(self):
        """Configure the argparser with arguments and defaults."""

        self._argparser.add_argument('-V','--version',
                action = 'version',
                version = '%(prog)s v{version}'.format(version=__version__))

        self._argparser.add_argument('-v','--verbose',
                help = """Be more verbose about each step.
                (default: %(default)s).""",
                action = 'store_true',
                dest = 'verbose',
                default = False,
                required = False)

        self._argparser.add_argument("-c", "--config-file",
                help="Specify PROST config file (default: %(default)s)",
                metavar = 'FILE',
                dest = 'config_file',
                default = DEFAULT_PROST_CONFIG_FILE,
                required = False)

        self._argparser.add_argument('--samples-filelist',
                help = """Name of the file containing the list of sample names
                          and filenames (default: %(default)s).""",
                metavar = 'FILE',
                action = 'store',
                dest = 'samples_filelist',
                default = DEFAULT_SAMPLES_FILELIST,
                required = False)

        self._argparser.add_argument('-d', '--output-dir',
                help=argparse.SUPPRESS,
                # help = """REWRITE - NOT IMPLEMENTED YET Name of directory will
                # be appended with date of run. Name of directory to store Prost!
                # output (default: %(default)s)""",
                metavar = 'DIR',
                action = 'store',
                dest = 'output_dir',
                default = DEFAULT_PROST_OUTPUT_DIR,
                required = False)

        self._argparser.add_argument('-s','--skip-sequence-alignments',
                help = """Skip all sequence alignments (e.g. BBMap searches
                against the genome and the annotation sequence databases)
                (default: %(default)s).""",
                action = 'store_true',
                dest = 'skip_sequence_alignments',
                default = False,
                required = False)

        self._argparser.add_argument('-sg','--skip-genome-alignment',
                help = """Skip the GenomeAlignment (e.g. BBMap search against the
                genome) only (default: %(default)s).""",
                action = 'store_true',
                dest = 'skip_genome_alignment',
                default = False,
                required = False)

        self._argparser.add_argument('-sa','--skip-annotation-alignments',
                help = """Skip the AnnotationAlignments (e.g. BBMap search
                against the annotation databases) only (default:
                %(default)s).""",
                action = 'store_true',
                dest = 'skip_annotation_alignments',
                default = False,
                required = False)

        self._argparser.add_argument('--min-seq-count',
                help = """The minimum number of times a sequence must appear
                (in total) across all samples to be included in the data set
                (default: %(default)s).""",
                metavar = 'NUM',
                type = int,
                action = 'store',
                dest = 'min_seq_count',
                default = DEFAULT_MIN_SEQ_COUNT,
                required = False)

        self._argparser.add_argument('--min-seq-length',
                help = """Sample sequences with length less than min-seq-length
                are ignored by prost (default: %(default)s).""",
                metavar = 'NUM',
                type = int,
                action = 'store',
                dest = 'min_seq_length',
                default = DEFAULT_MIN_SEQ_LENGTH,
                required = False)

        self._argparser.add_argument('--max-seq-length',
                help = """Sample sequences with length greater than
                max-seq-length are ignored by prost (default: %(default)s).""",
                metavar = 'NUM',
                type = int,
                action = 'store',
                dest = 'max_seq_length',
                default = DEFAULT_MAX_SEQ_LENGTH,
                required = False)

        self._argparser.add_argument('-w', '--wiggle',
                help = """The amount of genetic wiggle (default: %(default)s).""",
                metavar = 'NUM',
                type = int,
                action = 'store',
                dest = 'wiggle',
                default = DEFAULT_WIGGLE,
                required = False)

        self._argparser.add_argument('--max-threads',
                help = """The maximum number of threads (processors) to use for
                various operations (currently just aligners) (%(default)s).""",
                metavar = 'NUM',
                type = int,
                action = 'store',
                dest = 'max_threads',
                default = DEFAULT_MAX_THREADS,
                required = False)

        self._argparser.add_argument('--max-memory',
                help = """The maximum amount of memory (in GB) to allocate to
                the aligner.  If not set, defaults to 85 percent of system memory.""",
                metavar = 'NUM',
                type = int,
                action = 'store',
                dest = 'max_memory',
                default = get_85perc_system_memory(),
                required = False)

        self._argparser.add_argument('--max-locations-to-report',
                help = """The maximum number of genomic locations to report per
                sequence (if exceeded, it is simply reported as TML (too many
                locations) (default: %(default)s).""",
                metavar = 'NUM',
                type = int,
                action = 'store',
                dest = 'max_locations_to_report',
                default = DEFAULT_MAX_LOCATIONS_TO_REPORT,
                required = False)

        self._argparser.add_argument('--max-locations-allowed',
                help = """The maximum number of genomic locations allowed per
                sequence.  If exceeded for a given sequence, only the number of
                genomic locations is preserved, while the actual locations are
                erased, and the sequence is not binned (i.e. it is only
                included in the isomiRs or no_genomic_hits tab).  Reported as
                "MLAE" (maximum locations allowed exceeded).  This may help
                reduce Prost! memory consumption in some data sets.
                (default: %(default)s).""",
                metavar = 'NUM',
                type = int,
                action = 'store',
                dest = 'max_locations_allowed',
                default = DEFAULT_MAX_LOCATIONS_ALLOWED,
                required = False)

        self._argparser.add_argument('-C','--cache',
                help=argparse.SUPPRESS,
                # help = 'Use the cache (both read and write) (default: %(default)s).',
                action = 'store_true',
                dest = 'cache',
                default = False,
                required = False)

        self._argparser.add_argument('-Cr','--cache-read',
                help=argparse.SUPPRESS,
                # help = 'NOT YET IMPLEMENTED Read from the cache (default: %(default)s)',
                action = 'store_true',
                dest = 'cache_read',
                default = False,
                required = False)

        self._argparser.add_argument('-Cw','--cache-write',
                help=argparse.SUPPRESS,
                # help = 'NOT YET IMPLEMENTED Write to the cache (default: %(default)s)',
                action = 'store_true',
                dest = 'cache_write',
                default = False,
                required = False)

        self._argparser.add_argument('--output-file-prefix',
                help = """The prefix to be used for output filenames. (default:
                    %(default)s).""",
                metavar = 'FILE',
                action = 'store',
                dest = 'output_file_prefix',
                default = DEFAULT_OUTPUT_FILE_PREFIX,
                required = False)

        self._argparser.add_argument('--alignment-search-input-filename',
                help=argparse.SUPPRESS,
                # help = """Name of the generated search file used for the
                # Genome and Annotation Alignments (default: %(default)s).""",
                metavar = 'FILE',
                action = 'store',
                dest = 'input_file_for_alignments',
                default = DEFAULT_SEARCH_INPUT_FILE_FOR_ALIGNMENTS,
                required = False)

        self._argparser.add_argument('--mature-mir-annotation-fasta',
                help = """Optionally provide the mature miR annotation file and
                an extra column "MainSeqMatchesAnnotationFile" is added to the
                "by_annotation" tab.  See prost.config.example for more info.
                (default: %(default)s).""",
                metavar = 'FILE',
                action = 'store',
                dest = 'mature_mir_annotation_fasta',
                default = DEFAULT_MATURE_MIR_ANNOTATION_FASTA,
                required = False)

        self._argparser.add_argument('--create-tables',
                help=argparse.SUPPRESS,
                # help = 'Create the database caching tables.',
                action = 'store_true',
                dest = 'create_tables',
                default = False,
                required = False)

        self._argparser.add_argument('-Fr','--read-from-file-cache',
                help=argparse.SUPPRESS,
                # help = """If enabled, will try to load preexisting file
                # cache(s) from one or more stages (this is essentially a file
                # cache approach, in contrast to the database cache approach)
                # (default: %(default)s).""",
                action = 'store_true',
                dest = 'read_from_file_cache',
                default = False,
                required = False)

        self._argparser.add_argument('-Fw','--write-to-file-cache',
                help=argparse.SUPPRESS,
                # help = """If enabled, will write out to the file cache.
                # (default: %(default)s).""",
                action = 'store_true',
                dest = 'write_to_file_cache',
                default = False,
                required = False)

    def configure(self, db_proxy=None):
        # Set the argparser
        self._setup_argparser()

        # First parse, simply to get conf_file, in case it is not named
        # prost.config.
        args = self._argparser.parse_args()

        # Confirm existence of config file
        if not os.path.isfile(os.path.expanduser(args.config_file)):
            raise ConfigurationException, (
                "Configuration file {} does not exist.".format(
                    os.path.expanduser(args.config_file)))

        # Parse config file
        self._config.read(os.path.expanduser(args.config_file))

        # Combine defaults hash with GENERAL section of config file
        config_options = self.get_section('General')
        default_and_conffile_options = dict(
                list(self._defaults.items()) + list(config_options.items()))

        # Load defaults and config file options to argparser
        self._argparser.set_defaults(**default_and_conffile_options)

        # Second parse, to grab everything else.
        self.general = self._argparser.parse_args()

        # Lowercase the species for consistency
        self.general.species = self.general.species.lower()

        # Check the validity of the config file / command line options
        self._check_validity()

        # Create the output filenames
        self._create_output_filenames()

        ## Database Caching Stuff

        # Connect to optional database cache
        if db_proxy:
            self._connect_to_optional_database_cache(db_proxy)

        # Optionally create database tables
        self._create_database_tables()

        # If caching, create the genomic alignment seq_db entry in database if
        # needed
        self._create_genomic_alignment_seq_db_entry()

    def get_section(self, section_name):
        """Parse the passed config file section, converting particular fields
        to their known types (e.g. Boolean, Int, etc).

        Arguments:
            section_name: The name of a given config file section.
        Returns:
            A dict corresponding to configured options in the passed section.
        """

        conf = {}
        for field, val in self._config.items(section_name):
            if val == '':
                raise ConfigurationException, \
                    ("In prost.config, field '{}' is incorrectly blank.".format(
                        field))
            if field in CONFIG_FILE_BOOLEAN_FIELDS:
                conf[field] = self._config.getboolean(section_name, field)
            elif field in CONFIG_FILE_INT_FIELDS:
                conf[field] = self._config.getint(section_name, field)
            else:
                conf[field] = self._config.get(section_name, field)
        return conf

    def _create_database_tables(self):
        """Optionally create database tables."""
        if self.general.create_tables:
            if not self._config.has_section('Cache'):
                print ("WARNING: Ignoring request to create database tables "
                        "when no 'Cache' section has been specified.")
            else:
                SequenceDB.create_table()
                ShortSeqCache.create_table()
                GenomicLocationCache.create_table()

    def _create_genomic_alignment_seq_db_entry(self):
        """If caching, create (if it does not exist) the genomic alignment
        sequence database entry in the sequence_dbs database table.

        Note that this serves as a nice way to test the functioning of the db
        connection in the Configuration stage as well.
        """
        if self.general.cache:
            # Get or create the sequence database.
            conf_genome_alignment = self.get_section('GenomeAlignment')
            seq_db = conf_genome_alignment['db']
            try:
                SequenceDB.get(SequenceDB.name == seq_db)
            except DoesNotExist:
                SequenceDB.create(name = seq_db)

    def _check_validity(self):
        """Confirm the validity of the config file."""
        self._assert_general_section_has_required_options()
        self._assert_alignment_sections_have_required_options()
        self._assert_annotation_alignment_sections_have_additional_required_options()
        self._assert_alignment_sections_have_unique_names()
        self._assert_alignment_tools_are_supported()
        self._assert_annotation_alignment_types_are_supported()
        self._assert_optional_database_cache_is_valid()
        self._assert_non_negative_fields()
        self._assert_positive_fields()
        self._assert_sensible_configuration_values()

    def _create_output_filenames(self):
        """Create the output filenames from the user-specified file prefix and
        the hard-coded file suffixes."""

        if self.general.output_file_prefix:
            prefix = "{}_".format(self.general.output_file_prefix)
        else:
            prefix = ""

        self.general.output_file_isomirs = "{}{}.tsv".format(prefix,
                OUTPUT_FILE_SUFFIX_ISOMIRS)
        self.general.output_file_no_hits = "{}{}.tsv".format(prefix,
                OUTPUT_FILE_SUFFIX_NO_HITS)
        self.general.output_file_comp_by_gen_loc = "{}{}.tsv".format(prefix,
                OUTPUT_FILE_SUFFIX_COMP_BY_GEN_LOC)
        self.general.output_file_comp_by_seed = "{}{}.tsv".format(prefix,
                OUTPUT_FILE_SUFFIX_COMP_BY_SEED)
        self.general.output_file_comp_by_annotation = "{}{}.tsv".format(prefix,
                OUTPUT_FILE_SUFFIX_COMP_BY_ANNOTATION)
        self.general.output_file_mirror_mirs = "{}{}.tsv".format(prefix,
                OUTPUT_FILE_SUFFIX_MIRROR_MIRS)
        self.general.output_file_arm_switch = "{}{}.tsv".format(prefix,
                OUTPUT_FILE_SUFFIX_ARM_SWITCH)

        # Save all of the output filenames for later...
        self.general.output_files = [
            self.general.output_file_isomirs,
            self.general.output_file_comp_by_gen_loc,
            self.general.output_file_comp_by_annotation,
            self.general.output_file_comp_by_seed,
            self.general.output_file_mirror_mirs,
            self.general.output_file_arm_switch,
            self.general.output_file_no_hits,
        ]

    def _connect_to_optional_database_cache(self, db_proxy):
        """Setup the optional database connection."""
        try:
            conf_cache = self.get_section('Cache')

            # Set the database type
            db = conf_cache['type'].lower()
            if db == 'mysql':
                db_cls = MySQLDatabase
            elif db == 'postgres':
                db_cls = PostgresqlDatabase
            elif db == 'sqlite':
                db_cls = SqliteDatabase
            else:
                raise ControlFlowException, \
                        "ERR310: Shouldn't be possible to reach here."

            # Initialize the database
            if db_cls == SqliteDatabase:
                database = db_cls(conf_cache['database'])
            else:
                database = db_cls(conf_cache['database'],
                        user=conf_cache['username'],
                        passwd=conf_cache['password'])
            db_proxy.initialize(database)

        except ConfigParser.NoSectionError:
            """Section 'Cache' was not specified by user.  Skip."""
            pass

    def _assert_general_section_has_required_options(self):
        """The options in the 'General' section for the most part can be
        handled either via the command line or the config file.  There are
        a few options that we reserve just for the config file however.  Here
        they are, and the reasoning behind them:

        Required options in the [General] section:
            species: We require this in the config file simply because so much
                depends upon this field (especially in the annotation section).
        """

        for field in CONFIG_FILE_REQUIRED_GENERAL_FIELDS:
            if not self._config.has_option('General', field):
                raise ConfigurationException, \
                        ("Required field '{}' not found in the 'General' "
                         "section.".format(field))

    def _assert_alignment_sections_have_required_options(self):
        """All config file sections named '*Alignment*' must have required
        fields specified in CONFIG_FILE_REQUIRED_ALIGNMENT_FIELDS."""

        for sec_name, alignment in self.alignments.iteritems():
            for field in CONFIG_FILE_REQUIRED_ALIGNMENT_FIELDS:
                if field not in alignment.iterkeys():
                    raise ConfigurationException, \
                            ("Required field '{}' not found in alignment "
                             "section '{}'".format(field, sec_name))

    def _assert_annotation_alignment_sections_have_additional_required_options(self):
        """All config file sections named "AnnotationAlignment*" must have
        additional required fields as specified in
        CONFIG_FILE_ADDITIONAL_REQUIRED_ANNOTATION_ALIGNMENT_FIELDS.  These are
        in additionl to those found in
        CONFIG_FILE_REQUIRED_ALIGNMENT_FIELDS."""

        for sec_name, alignment in self.alignments.iteritems():
            if not re.search(r'^AnnotationAlignment', sec_name):
                continue
            for field in CONFIG_FILE_ADDITIONAL_REQUIRED_ANNOTATION_ALIGNMENT_FIELDS:
                if field not in alignment.iterkeys():
                    raise ConfigurationException, \
                            ("Required field '{}' not found in annotation "
                             "alignment section '{}'".format(field, sec_name))

    def _assert_alignment_sections_have_unique_names(self):
        """The 'name' field must be unique among sections name '*Alignment*'."""

        names = {}
        for sec_name, alignment in self.alignments.iteritems():
            if alignment['name'] in names.itervalues():
                raise ConfigurationException, \
                        ("Alignment names must be unique (name '{}' in section "
                         "'{}' is not unique.".
                         format(alignment['name'], sec_name))
            else:
                names[sec_name] = alignment['name']

    def _assert_alignment_tools_are_supported(self):
        """Confirm alignment tools are supported."""

        for sec_name, alignment in self.alignments.iteritems():
            if alignment['tool'].lower() not in CONFIG_FILE_SUPPORTED_ALIGNMENT_TOOLS:
                raise ConfigurationException, \
                        ("Alignment tool '{}' in section '{}' not supported.\n"
                         "We support these alignment tools: {}.".
                         format(alignment['tool'], sec_name,
                             ", ".join(CONFIG_FILE_SUPPORTED_ALIGNMENT_TOOLS)))

    def _assert_annotation_alignment_types_are_supported(self):
        """Confirm the annotation alignment types (e.g
        MirbaseHairpinAnnotation) are supported."""

        for sec_name, alignment in self.alignments.iteritems():
            if not re.search(r'^AnnotationAlignment', sec_name):
                continue
            if alignment['type'] not in \
                    CONFIG_FILE_SUPPORTED_ANNOTATION_ALIGNMENT_TYPES:
                raise ConfigurationException, \
                        ("Annotation Alignment type '{}' in section '{}' is "
                         "not supported.\nWe support these types: {}.".
                         format(alignment['tool'], sec_name,
                             ", ".join(
                                 CONFIG_FILE_SUPPORTED_ANNOTATION_ALIGNMENT_TYPES
                                 )))

    def _assert_optional_database_cache_is_valid(self):
        """Confirm the 'Cache' section is valid."""

        try:
            cache_config = dict(self._config.items('Cache'))

            # Are the required fields present in the 'Cache' section?
            for field in CONFIG_FILE_REQUIRED_CACHE_FIELDS:
                if field not in cache_config.iterkeys():
                    raise ConfigurationException, \
                            ("Required field '{}' not found in 'Cache' section.".
                             format(field))

            # Is the database supported?
            if cache_config['type'] not in CONFIG_FILE_SUPPORTED_DATABASES:
                raise ConfigurationException, \
                        ("Database '{}' in section 'Cache' is not supported.\n"
                         "We support these databases: {}".
                         format(cache_config['type'],
                             ", ".join(CONFIG_FILE_SUPPORTED_DATABASES)))

        except ConfigParser.NoSectionError:
            """Section 'Cache' was not specified by user.  Skip."""
            pass

    def _assert_non_negative_fields(self):
        """Ensure that certain fields are indeed non-negative.

        Note that this is applied to the merged product of the config file and
        the command line parameters.
        """

        for section_name in self._config.sections():
            for field in CONFIGURATION_NON_NEGATIVE_FIELDS:
                # Set val to arbitrary large positive number
                val = ONE_MILLION

                if section_name == 'General':
                    try:
                        val = getattr(self.general, field)
                    except AttributeError:
                        pass
                else:
                    # 'field' is not in General section, try Alignment secs.
                    try:
                        val = self._config.get(section_name, field)
                    except ConfigParser.NoOptionError:
                        pass
                if val < 0:
                    if section_name != "General":
                        in_section_str = " (in section '{}')".format(section_name)
                    else:
                        in_section_str = ""
                    raise ConfigurationException, \
                        ("Field {}{} cannot be negative (currently {}).".format(
                            field, in_section_str, val))

    def _assert_positive_fields(self):
        """Ensure that certain fields are indeed positive.

        Note that this is applied to the merged product of the config file and
        the command line parameters.
        """

        for section_name in self._config.sections():
            for field in CONFIGURATION_POSITIVE_FIELDS:
                # Set val to arbitrary large positive number
                val = ONE_MILLION

                if section_name == 'General':
                    try:
                        val = getattr(self.general, field)
                    except AttributeError:
                        pass
                else:
                    # 'field' is not in General section, try Alignment secs.
                    try:
                        val = self._config.get(section_name, field)
                    except ConfigParser.NoOptionError:
                        pass
                if val <= 0:
                    if section_name != "General":
                        in_section_str = " (in section '{}')".format(section_name)
                    else:
                        in_section_str = ""
                    raise ConfigurationException, \
                        ("Field {}{} must be positive (currently {}).".format(
                            field, in_section_str, val))

    def _assert_sensible_configuration_values(self):
        """Ensure that the user has chosen sensible configuration values."""

        # 1. min_seq_len must be less than max_seq_len
        if self.general.max_seq_length <= self.general.min_seq_length:
            raise ConfigurationException, \
                ("min_seq_length must be less than max_seq_length.")

        # 2. max_locations_allowed must be >= max_locations_to_report
        if (self.general.max_locations_to_report >
                    self.general.max_locations_allowed):
            raise ConfigurationException, \
                ("max_locations_to_report cannot be greater than "
                 "max_locations_allowed.")


###############
### Classes ###
###############

class Samples(collections.OrderedDict):
    """Represents the list of samples.

    key: Name of the sample, e.g. 'sampleA'.
    val: Name of the fasta file holding the sample data."""

    def _read_samples_filelist(self, filelist):
        """Read in the samples file list and populate self's dict."""
        with open(os.path.expanduser(filelist), 'rU') as f:
            for sample_line in f:
                if sample_line[0] == '#':
                    continue
                name_and_file = sample_line.strip().split()
                name_and_file.reverse()
                if (len(name_and_file) != 2):
                    raise ConfigurationException, \
                        ("The file '{}' requires exactly two entries per line "
                         "(sample filename and sample name)".format(filelist))
                [sample_name, sample_file] = name_and_file
                if sample_name in self.iterkeys():
                    raise ConfigurationException, \
                        ("The sample name '{}' is illegally "
                         "duplicated in '{}'.".format(sample_name, filelist))
                self[sample_name] = sample_file


class ShortSeq(SlotPickleMixin):
    """Represents a single short RNA sequence.

    Attributes:
        seq_str (str): The string (nucleotides) representation of this ShortSeq.
        designation_integer (int): The integer part of the designation of this
            ShortSeq (i.e. left of the decimal).
        sum_norm (float): The sum of the normalized read counts for this
            ShortSeq across all samples.
        ambiguous (bool): Sequences are defined as 'ambiguous' if they are
            within_wiggle of more than one GenLocBin bin_starter.
        cached (bool): Whether or not this ShortSeq has been cached in the DB.
        count_genomic_locations (int): The number of genomic locations.
            Recorded in case it exceeds max_locations_allowed (in that case,
            genomic_locations[] is emptied, so this field is our record of the
            actual number of locations.)
        count_no_hit_genomic_locations (int): The number of "no_hit" genomic
            locations (see Arg 'no_hit_genomic_locations' below). Recorded in
            case it exceeds max_locations_allowed (in that case,
            no_hit_genomic_locations[] is emptied, so this field is our record
            of the actual number of locations.)
        max_locations_allowed_exceeded (bool): True it the number of genomic
            locations for this ShortSeq exceeds max_locations_allowed.
        samples_counts (SampleCountsWithNorms): An object representing the
            read counts for this ShortSeq across samples.  Also holds the
            normalized read counts for this ShortSeq across samples.
        wiggle_fingerprints ((int, (int,), (int,)), (int, (int,), (int,))): The
            two wiggle fingerprints of this ShortSeq. See
            to_wiggle_fingerprints().
        genomic_locations ([AlignmentExecutionHit]): A list of genomic locations
            at which this ShortSeq is found.  We use the AlignmentExecutionHit
            as the representation of a short sequence's Genomic Location; in
            this context, the AlignmentExecutionHit's reference name, start, and
            stop represent the Genomic Location of the short sequence.
        no_hit_genomic_locations ([AlignmentExecutionHit]): A list of genomic
            locations at which this ShortSeq is found but for which Prost has
            post-filtered out (e.g. for having too many 3p mismatches).
        annotations ([Annotations]): The list of annotations for this ShortSeq.

    """
    __slots__ = (
            'seq_str',
            'designation_integer',
            '_sum_norm',
            'ambiguous',
            'cached',
            'count_genomic_locations',
            'count_no_hit_genomic_locations',
            'max_locations_allowed_exceeded',
            'samples_counts',
            'genomic_locations',
            'no_hit_genomic_locations',
            'wiggle_fingerprints',
            'annotations',
            'iso_5p',
            'iso_3p',
            'iso_add',
            'iso_snp_seed',
            'iso_snp_central_offset',
            'iso_snp_central',
            'iso_snp_supp',
            'iso_snp',
            #'memo_seed_shifted',
            #'memo_seed_edited',
            #'memo_3p_supplementary_edited',
            #'memo_3p_modifications',
            #'memo_3p_untemplated_addition',
            #'memo_5p_untemplated_or_edited',
            #'memo_other_edited',
    )

    @property
    def sum_norm(self):
        """The sum of the normalized read counts for this ShortSeq."""
        return self.samples_counts.sum_norm()

    @property
    def is_one(self):
        """Is this short sequence designated a 'one'?"""
        return (self.designation_integer == DESIGNATION_ONE)

    @property
    def is_no_hit(self):
        """Did this ShortSeq have no hits against the genome?"""
        return (self.designation_integer == DESIGNATION_NO_HIT)

    @property
    def has_indelnts(self):
        """Does this ShortSeq have any genomic locations with gaps?"""
        for gl in self.genomic_locations:
            if gl.designation:
                if re.search(r'\.{}$'.format(INDEL_INDICATOR), gl.designation):
                    return True
        return False

    @property
    def seed(self):
        """The seed region of this ShortSeq.

        Returns:
            str: The 7 nucleotide seed region (nts 2-8).

        """
        return self.seq_str[1:8]

    def __init__(self, seq_str, num_samples):
        self.seq_str = seq_str
        self._sum_norm = -1.0
        self.cached = False
        self.designation_integer = None
        self.ambiguous = None
        self.count_genomic_locations = None
        self.count_no_hit_genomic_locations = None
        self.max_locations_allowed_exceeded = False
        self.wiggle_fingerprints = None
        self.samples_counts = SamplesCountsWithNorms(num_samples)
        self.genomic_locations = []
        self.no_hit_genomic_locations = []
        self.annotations = []
        self.iso_5p = 'NA'
        self.iso_3p = 'NA'
        self.iso_add = 'NA'
        self.iso_snp_seed = False
        self.iso_snp_central_offset = False
        self.iso_snp_central = False
        self.iso_snp_supp = False
        self.iso_snp = False

    def __repr__(self):
        """What to return when pretty printing."""
        #if self.genomic_locations:
        #    gen_locs = self.genomic_locations[0:2]
        #else:
        #    gen_locs = None
        return pprint.pformat({
            "seq": self.seq_str,
            # "cached": self.cached,
            # "designation_integer": self.designation_integer,
            # "sum_norm": self.sum_norm,
            # "(first two) gen_locs": gen_locs,
            # "annotations": self.annotations
            })

    def __str__(self):
        """What to return when printing."""
        return self.__repr__()

    def set_ambiguous(self, ambiguous_to):
        """Set the list of ambiguous sequences for this ShortSeq.  A ShortSeq
        is ambiguous if it matches up to more than one 'one_bin_starter' or
        more than on 'non_one_bin_starter'. (Note that self.ambiguous is a list
        of sequence strings, *not* ShortSeq objects.)

        The list is sorted alphabetically.
        """
        assert isinstance(ambiguous_to, list)
        self.ambiguous = sorted(ambiguous_to)

    def set_wiggle_fingerprints(self, max_seq_length):
        """Set the 'Wiggle Fingerprints' (wiggle_fingerprints) for this short
        sequence.

        See to_wiggle_fingerprints() for more information.

        wiggle_fingerprints are used to speed within_wiggle searches.

        Returns:
            ((int, (int,), (int,)), (int, (int,), (int,))): Two wiggle
                fingerprints. See to_wiggle_fingerprints().
        """
        assert (len(self.genomic_locations) != 0 or self.is_no_hit), \
                "ShortSeq.genomic_locations[] cannot be empty when executing "\
                "set_wiggle_fingerprints() unless it is a no hit designation."
        self.wiggle_fingerprints = to_wiggle_fingerprints(
            self.genomic_locations, max_seq_length)
        return self.wiggle_fingerprints

    def within_wiggle(self, other, wiggle):
        """Determines if this ShortSeq is within genetic wiggle of the 'other'
        ShortSeq.

        Genetic Wiggle is defined as follows.

        Sequence A can be at multiple genomic locations.  Sequence B can be at
        multiple genomic locations as well.  Sequence B will be thrown in the same
        bin as Sequence A if it is within the "genetic wiggle" of A.

        For Sequence B to be within the wiggle, it's basically got to overlap on
        exactly the same genetic locations (no more, no less) as Sequence A, with
        the A.start +/- 'x' bps from B.start, and A.stop +/- 'x' bps from B.stop.

        Another way for Sequence B to be within wiggle of Sequence A is for
        Sequence B to be entirely contained within Sequence B at all genetic
        locations.

        We define a ShortSeq to *not* be within wiggle of itself.

        Assumption: minimum length of any sequence is >= (2*wiggle + 1).

        WARNING: This method is NOT symmetric.  A may be within_wiggle of B,
        but B is then not nececssarily within_wiggle of A.

        Arguments:
            other: The other ShortSeq.
            wiggle: The size of the wiggle.
        """

        within = True

        if self is other:
            # A ShortSeq is not within wiggle of itself.
            within = False
        elif(len(self.genomic_locations) != len(other.genomic_locations)):
            # Cannot have different number of genomic locations and be within
            # wiggle.
            within = False
        else:
            # Ensure each pair of genomic locations is within wiggle.
            gen_loc_pairs = zip(self.genomic_locations, other.genomic_locations)
            for gen_loc_self, gen_loc_other in gen_loc_pairs:
                if not gen_loc_self.within_wiggle(gen_loc_other, wiggle):
                    within = False
                    break

        return within


class SamplesCounts(SlotPickleMixin):
    """A single SamplesCounts object generically tracks some per-sample count
    for each Sample in Samples.

    A single SamplesCounts object may be attached to a ShortSeq, a Bin, or a
    Bins object.  For a given ShortSeq object, the associated SamplesCounts
    object tracks its read counts (absolute).  For a Bin object, the associated
    SamplesCounts object tracks the per-sample read totals for all the Bin's
    ShortSeq objects.  Bins has a similar relationship to Bin in this regard.
    """

    __slots__ = ('num_samples', '_counts',)

    def __init__(self, num_samples):
        self.num_samples = num_samples
        # Arrays for more efficient storage.
        self._counts = array('L', (0 for i in xrange(num_samples)))

    def __iadd__(self, other):
        """Operator overload of '+='.

        Args:
            other (SamplesCounts): The other SampleSeqCounts to which you
                wish to add to this SampleSeqCounts.
        """
        for i in xrange(self.num_samples):
            self._counts[i] += other._counts[i]
        return self

    def __repr__(self):
        return pprint.pformat(self._counts)

    def incremenent_count(self, sample_num):
        """Increment the 'sample_num'th sample count (absolute)."""
        self._counts[sample_num] += 1

    def total_count(self):
        """Calculates the total absolute count across samples."""
        return sum(self._counts)

    def count(self, sample_num):
        """Getter for the (absolute) count of a given sample (i.e.
        sample_num)."""
        return self._counts[sample_num]

    def counts(self):
        """Generator of the (absolute) counts for all samples."""
        for i in xrange(len(self._counts)):
            yield self._counts[i]

    def na(self):
        """Return an array of 'NA's for outputting purposes."""
        return ['NA'] * self.num_samples

    def perc(self, samples_totals):
        """Calculate per-sample percents, where 'self' is the numerator and
        'samples_totals' is the denominator.  Used for tsv output.

        Args:
            samples_totals (SamplesCounts): The totals with which to divide
                'self' by.

        Returns:
            [floats]: An array of floats representing per-sample
                percentages.
        """
        percents = []
        for i in xrange(self.num_samples):
            if samples_totals._counts[i] != 0:
                percents.append(float(self._counts[i]) / samples_totals._counts[i])
            else:
                percents.append(0.0)
        return percents


class SamplesCountsWithNorms(SamplesCounts):
    """Extension of SamplesCounts which also caches per-sample normalizations.

    The normalized counts should be computed once all the relevant attached
    objects (e.g. all ShortSeq objects) have had their SamplesCounts objects
    attached and absolute counts tallied.
    """

    __slots__ = ('_norms',)

    def __init__(self, num_samples):
        super(SamplesCountsWithNorms, self).__init__(num_samples)
        # Arrays for more efficient storage.

        self._norms = array('d', (0 for i in xrange(num_samples)))

    def __repr__(self):
        return pprint.pformat((self._counts, self._norms))

    def sum_norm(self):
        """Returns the sum of the normalized read counts across all samples for
        this SamplesCounts."""
        return sum(self._norms)

    def normalize_sample(self, sample_num, sample_total):
        """Normalize the sample 'sample_num' with the total sequence count
        across the sample (i.e. 'sample_total')."""
        if sample_total == 0:
            self._norms[sample_num] = 0.0
        else:
            self._norms[sample_num] = \
                    self._counts[sample_num] * ONE_MILLION / sample_total

    def normalize_samples(self, samples_totals):
        """Normalize all samples using samples_totals, which is a SamplesCounts
        object containing the per-sample totals.

        Args:
            samples_totals (SamplesCounts): The per-sample totals.

        """
        for i in xrange(self.num_samples):
            self.normalize_sample(i, samples_totals._counts[i])

    def norm(self, sample_num):
        """Getter for the normalized count of a given sample (i.e.
        sample_num)."""
        return self._norms[sample_num]

    def norms(self):
        """Generator of the normalized counts for all samples."""
        for i in xrange(len(self._norms)):
            yield self._norms[i]


class ShortSeqs(dict):
    """Container for short RNA sequences.

    This is a singleton object container to hold unique short RNA sequences (of
    type ShortSeq).  It functions as a dictionary, where each key is the sequence
    string of a ShortSeq object and each value is the corresponding ShortSeq
    object.  This singleton is used extensively throughout the code base.
    """

    def __init__(self, configuration, *args, **kwargs):
        super(ShortSeqs, self).__init__(self, *args, **kwargs)
        self._conf = self.setup_configuration(configuration)
        self._designation_index = {}
        self._wiggle_fingerprint_bins_low = {}
        self._wiggle_fingerprint_bins_high = {}
        self.per_sample_read_totals = None

        # Memo to determine when ShortSeqs' first ShortSeq has been designated.
        self._first_sequence_has_been_designated = False

    def setup_configuration(self, config):
        """Store only the configuration items that are needed by ShortSeqs.

        This is done instead of storing a reference to the Configuration
        singleton because the Configuration singleton contains a reference to
        the argparse object, which cannot be pickled.  (The database reference
        contained in the Configuration singleton probably also cannot be
        pickled.)

        Arguments:
            config - The Configuration singleton.
        """
        conf = {}
        conf['species'] = config.general.species
        conf['min_seq_length'] = config.general.min_seq_length
        conf['max_seq_length'] = config.general.max_seq_length
        conf['min_seq_count'] = config.general.min_seq_count
        conf['wiggle'] = config.general.wiggle
        return conf

    def process_sample_fastas(self, samples):
        """Process the fasta files containing the samples.

        Specifically, this method reads the fasta files, normalizes them per
        sample file, and rejects sequences with low reads.

        Args:
            samples (Samples): The (populated) Samples singleton object.
        """

        self.read_sample_fastas(samples)
        self.reject_low_read_seqs()

    def read_sample_fastas(self, samples):
        """Read FASTA file of samples and create dictionary with associated
        (absolute only for now) counts per sample.

        Args:
            samples (Samples): The (populated) Samples singleton object.
        """

        pmsg("Reading in all fasta files...")
        num_samples = len(samples)

        for sample_num, sample_name in enumerate(samples):
            num_entries = sum(1 for line in open(samples[sample_name], 'rU'))/2
            filename = os.path.basename(samples[sample_name])
            progress = Progress(
                    'reading file{} ({})'.format(sample_num+1, filename),
                    100000, num_entries)

            # populate absolute counts for this sample
            min_seq_length = self._conf['min_seq_length']
            max_seq_length = self._conf['max_seq_length']
            with open(samples[sample_name], 'rU') as in_fasta:
                for line in in_fasta:
                    if line[0] == '>':
                        continue
                    progress.progress()
                    rec = line.rstrip()
                    if ((len(rec) >= min_seq_length) and
                            (len(rec) <= max_seq_length)):
                        try:
                            short_seq = self[rec]
                        except KeyError:
                            short_seq = self[rec] = ShortSeq(rec, num_samples)
                        short_seq.samples_counts.incremenent_count(sample_num)

            progress.done()

    def reject_low_read_seqs(self):
        """Reject sequences with low reads.

        Currently, if the sum of the reads across samples for a given sequence
        is less than --min-seq-count, then the sequence is rejected.

        In the future, we may automate this so that --min-seq-count is
        chosen based upon the character of the input fasta files."""

        progress = Progress("Rejecting low reads seqs")
        for seq_str in self.keys():
            short_seq = self[seq_str]
            if (short_seq.samples_counts.total_count()
                    < self._conf['min_seq_count']):
                del(self[seq_str])
        progress.done()

    def generate_fasta_search_file(self, cache, input_file_for_alignments):
        """Generate the fasta search file which will be used for alignments.
        """
        progress = Progress("Writing fasta search file")

        if cache:
            seq_strs_for_alignments = (short_seq.seq_str for short_seq in
                    self.itervalues() if short_seq.cached == False)
        else:
            seq_strs_for_alignments = self.iterkeys()

        with open(input_file_for_alignments, 'w') as fasta_search:
            for seq_str in seq_strs_for_alignments:
                fasta_search.write(">{}\n{}\n".format(seq_str, seq_str))

        progress.done()


    def designate_sequence(self, seq_str, designation_integer,
            max_locations_allowed=None, gen_locs=None, no_hit_gen_locs=None):
        """Bookkeeping to add the designation and genomic locations of short
        sequence 'seq_str' to the ShortSeqs and appropriate ShortSeq
        objects.

        The designtation is tracked by self._designation_index[], while the
        genomic locations is tracked within the relevant ShortSeq object
        contained within self's dict.  Additionally, we track the "Wiggle
        Fingerprints" for quicker 'within_wiggle*' calculations later.

        Finally, we also set the sequence's designtation_integer.

        This method should be called only once for each sequence.

        Args:
            seq_str (str): The short sequence to be designated and updated.
            designation_integer (int): The integer part of the designation
                of this short sequence (i.e. left of the decimal).
            max_locations_allowed (int): Optional - The maximum number of
                genomic locations allowed per sequence.  If exceeded, all
                genomic locations are erased for that sequence (count is still
                preserved).
            gen_locs ([GenomicLocation]): Optional - Locations on genome that
                both a) aligned to the sequence and b) have the designation
                'designation'.  This criteria means that some alignment hits
                are thrown away. This optional argument should only be included
                once per ShortSeq.
            no_hit_gen_locs ([GenomicLocation]): Optional - Locations on genome that
                aligned to the genome but did not pass Prost-post-filtering.
        """
        # setup
        short_seq = self[seq_str]
        assert isinstance(designation_integer, int)

        if no_hit_gen_locs:
            # Record number of no_hit_gen_locs:
            short_seq.count_no_hit_genomic_locations = len(no_hit_gen_locs)

            if len(no_hit_gen_locs) > max_locations_allowed:
                # Too many locations to save, toss them all to save memory.
                short_seq.no_hit_genomic_locations[:] = []
            else:
                # This is the only place that we sort no_hit genomic locations.
                short_seq.no_hit_genomic_locations[:] = sorted(no_hit_gen_locs)

        if gen_locs:
            # We are only setting gen_locs once for each sequence, right?
            assert len(short_seq.genomic_locations) == 0

            # Record number of gen_locs:
            short_seq.count_genomic_locations = len(gen_locs)

            if len(gen_locs) > max_locations_allowed:
                # Too many locations to save, toss them all to save memory.
                short_seq.genomic_locations[:] = []
                # Note it:
                short_seq.max_locations_allowed_exceeded = True
                # Reset designation appropriately
                designation_integer = \
                        DESIGNATION_MAX_LOCATIONS_ALLOWED_EXCEEDED
            else:
                # Reallocating list in some cases here is wee bit inefficient.
                # This is the only place that we sort no_hit genomic locations.
                short_seq.genomic_locations[:] = sorted(gen_locs)

                # update wiggle_fingerprint_bins
                (fp_low, fp_high) = \
                        short_seq.set_wiggle_fingerprints(
                                self._conf['max_seq_length'])

                self._wiggle_fingerprint_bins_low[fp_low] = \
                        self._wiggle_fingerprint_bins_low.get(fp_low, set())
                self._wiggle_fingerprint_bins_low[fp_low].add(short_seq)

                self._wiggle_fingerprint_bins_high[fp_high] = \
                        self._wiggle_fingerprint_bins_high.get(fp_high, set())
                self._wiggle_fingerprint_bins_high[fp_high].add(short_seq)

        # Set the designation_integer, which may be '12' (i.e. one or two).
        short_seq.designation_integer = designation_integer

        if (designation_integer != DESIGNATION_ONE_OR_TWO):
            # update designation_bins if we have our final designation.
            di = designation_integer
            self._designation_index[di] = self._designation_index.get(di, [])
            self._designation_index[di].append(seq_str)

    def first_sequence_designated_one_has_been_designated(self):
        """Memoized function to determine if the first sequence has been
        designated yet or not.
        """
        if self._first_sequence_has_been_designated:
            return True
        else:
            if (DESIGNATION_ONE in self._designation_index
                    and len(self._designation_index[DESIGNATION_ONE])):
                self._first_sequence_has_been_designated = True
                return True
            else:
                return False

    def ones(self):
        """Generator of all ones.

        Yields:
            ShortSeq: A ShortSeq object that has a DESIGNATION of 'one'.
        """
        for short_seq in self.itervalues():
            if short_seq.is_one:
                yield short_seq

    def non_ones_sorted_by_designation_then_sum_norm(self):
        """Returns a list of all sequences not designated as 'one', sorted
        first by designation (ascending) then by the sum of normalized counts
        (descending).

        Note:
            For reproducability, the list is finally sorted by the sequence
            string (descending), though this sorting should not affect the
            algorithm).

        Returns:
            [ShortSeq]: A list of ShortSeq objects not designated as 'one',
            sorted first by designation (ascending) then by maximum normalized
            count (descending), then by sequence string.
        """
        non_ones = (short_seq for short_seq in
                self.itervalues() if short_seq.designation_integer !=
                DESIGNATION_ONE)

        return sorted(non_ones, key=lambda (short_seq): (
            short_seq.designation_integer,
            -short_seq.sum_norm,
            short_seq.seq_str))

    def sorted_by_designation_then_sum_norm(self):
        """Returns a list of all sequences sorted first by designation
        (ascending) then by the sum of normalized counts (descending).

        Note: For reproducability, the list is finally sorted by the sequence
            string (descending), though this sorting should not affect the
            algorithm.

        Returns:
            [ShortSeq]: A sorted list of ShortSeq objects.
        """

        return sorted(self.itervalues(),
                key=lambda (short_seq): (
                    short_seq.designation_integer,
                    -short_seq.sum_norm,
                    short_seq.seq_str))

    def candidate_wigglers(self, short_seq):
        """Returns the set() of short seqs that are potentially within_wiggle
        to short_seq.

        Implementation detail: ...

        Good example as to why we need to check high/low weird combinations.
        Document!!!!
        #                     binlow   binhigh
        # starts1 = (58)   -> (40) ,    (60)
        # starts2 = (70)   -> (60) ,    (80)
        """

        # Get the two fingerprints
        (fp_low, fp_high) = short_seq.wiggle_fingerprints

        # The set of candidate short seqs is the union of the wiggle
        # fingerprint bins key'ed off the low and high fingerprints.
        return (self._wiggle_fingerprint_bins_low.get(fp_low, set())
                | self._wiggle_fingerprint_bins_high.get(fp_high, set())
                | self._wiggle_fingerprint_bins_high.get(fp_low, set())
                | self._wiggle_fingerprint_bins_low.get(fp_high, set()))

    def within_wiggle_of(self, short_seq, inverse):
        """Generator of ShortSeq objects that are within_wiggle of short_seq.

        Generator of all short sequences for which:

            a) each short sequence is within wiggle of the passed 'short_seq'
            b) the passed 'short_seq' is within wiggle of each short sequence

        When 'inverse' is True, we perform 'b' above; when 'inverse' is False,
        we perform 'a' above.

        We define a sequence to never be within_wiggle of itself.

        Args:
            short_seq (ShortSeq): A ShortSeq object.

        Yields:
            ShortSeq: A ShortSeq object.

        """
        if not short_seq.genomic_locations:
            # A no-hit ShortSeq cannot be within_wiggle of any others.
            return

        for other in self.candidate_wigglers(short_seq):
            if short_seq is other:
                # Skip identity.
                continue
            if not inverse:
                # We are generating seqs within_wiggle of short_seq
                if other.within_wiggle(short_seq, self._conf['wiggle']):
                    yield other
            else:
                # We are generating seqs for which short_seq is within_wiggle
                # of each of those seqs.
                if short_seq.within_wiggle(other, self._conf['wiggle']):
                    yield other

    def is_within_wiggle_of_ones(self, short_seq):
        """Is the passed ShortSeq within_wiggle of any ones?

        Args:
            short_seq (ShortSeq): A ShortSeq that we wish to test to determine
                if it is within_wiggle of any ones.
        Returns:
            bool: True if the short_seq is within_wiggle of at least one
                sequence designated 'one', False otherwise.
        """

        for other in self.candidate_wigglers(short_seq):
            if not other.is_one:
                continue
            if short_seq.within_wiggle(other, self._conf['wiggle']):
                return True
        return False

    def designation_step_one(self, alignment, max_locations_allowed):
        """Designates each ShortSeq object as either 3 or above, or "1_or_2",
        or "no_hit".

        This is the first pass of designation.  In this pass, no wiggle
        calculations are performed.  Instead, this pass merely looks as how well
        each sequence aligns to the genome.  All 3's and above as well as any
        sequence which does not hit to the genome can be definitively
        designated.

        In contrast, it is impossible to tell at this stage whether or not the a
        sequence is a one or a two; for such sequences, this method assigns a
        designation of DESIGNATION_ONE_OR_TWO.  See designation_step_two() for
        the stage in which ones and twos are definitively designated.

        Args:
            alignment (GenomeAlignment): A GenomeAlignment object.
            max_locations_allowed (int): Optional - The maximum number of
                genomic locations allowed per sequence.  If exceeded, all
                genomic locations are erased for that sequence (count is still
                preserved).

        """
        progress = Progress("Designation step ONE")
        logging.info("Designation step ONE...")
        sys.stderr.write('\n')

        assert isinstance(alignment, GenomeAlignment), \
                "Alignment {} is not a GenomeAlignment".format(alignment)

        # designate all hits (that have not been post-filtered by Prost)
        for hit_set in alignment.alignment_execution.single_sequence_hit_sets():
            if len(hit_set):
                query_sequence = hit_set[0].query_sequence
                try:
                    self[query_sequence]
                except KeyError:
                    # In cases where we are using alignment results that were
                    # generated outside of the current execution of Prost
                    # (perhaps with less stringent min seq count), we may find
                    # hits whose query sequences have been removed from
                    # ShortSeqs due to low seq count.  In that case, simply
                    # continue.
                    continue
                self.designate_genome_alignment_hit_set(hit_set,
                        alignment.alignment_execution.max_non_3p_mismatches,
                        alignment.alignment_execution.indelnt_penalty_multiplier,
                        max_locations_allowed)

        progress.done()

    def designate_genome_alignment_hit_set(self, hit_set, max_non_3p_mismatches,
            indelnt_penalty_multiplier, max_locations_allowed):
        """Designate all hits within a hit_set.

        A hit_set is a set of alignment hits of a single sequence against
        a reference genome.

        This method designates the hits within a hit_set, keeping only hits
        with the best sequence designtation.

        For example, while a hit_set may contain a mixture of hits with
        designtations one, two, and three, only the hits with a designation of one
        are kept.

        Note that only the integer part of the designation is considered (e.g.
        a designation '3.2' is considered equal to a designation '3.1').

        This method is the first pass. It designates 3's, 4's, 5's, etc (those
        will never change designation) and waffles on 1's and 2's, assigning
        those a designation of DESIGNATION_ONE_OR_TWO, meaning 'either a 1 or
        a 2'.

        Args:
            hit_set ([AlignmentExecutionHit]): A list of alignment hits which
                all contain the same query sequence.
            max_non_3p_mismatches (int): The maximum number of non 3p mismatches
                allowed per alignment.
            indelnt_penalty_multiplier (int): As it sounds.
            max_locations_allowed (int): Optional - The maximum number of
                genomic locations allowed per sequence.  If exceeded, all
                genomic locations are erased for that sequence (count is still
                preserved).

        """
        query_sequence = hit_set[0].query_sequence

        # Do not perform designation step one if already cached.
        if self[query_sequence].cached:
            return

        gen_locs = []
        no_hit_gen_locs = []
        lowest_designation = 100
        longest_hit_len = 0

        for hit in hit_set:
            # check if hit is a one
            curr_designation = 0

            if hit.is_no_hit:

                # This is a prost-post-filtered hit (e.g. too many 3p mismatches)
                # Store for now for potential inclusion in no_gen_hits excel tab.
                no_hit_gen_locs.append(hit)
            elif(hit.has_100_percent_core_identity and not hit.has_5p_mismatch):

                # Designation: 1, 2, or 3. (Is 100% ident.)
                if(hit.is_full_length):
                    # Designation: 1 or 2. (Is 100% ident and full length.)
                    gen_locs.append(hit)
                    lowest_designation = 2
                else:
                    # Designation: 3. (Is 100% ident but not full length.)
                    curr_designation = 3
                    if(lowest_designation < curr_designation):
                        # Only take the best designations.
                        continue
                    elif(hit.alignment_length < longest_hit_len):
                        # Only take the longest hits.
                        continue
                    elif(lowest_designation > curr_designation):
                        # Throw away any previous hits which have a designation
                        # larger than curr_designation (i.e. '3' in this case).
                        gen_locs[:] = []

                    # We accept this hit (for now).
                    longest_hit_len = hit.alignment_length

                    hit.designation = hit.calculate_designation_full(DESIGNATION_THREE)
                    gen_locs.append(hit)
                    lowest_designation = curr_designation

            else:

                # Designation: >= 4. (Not 100% ident, i.e. not a perfect hit.
                # May or may not be full length.)
                curr_designation = (
                        3
                        + hit.num_non_3p_mismatches
                        + (hit.num_indelnts * indelnt_penalty_multiplier)
                )
                # Set the maximum possible designation
                max_designation = 3 + max_non_3p_mismatches

                if(lowest_designation < curr_designation):
                    # Only take the best designations.
                    continue
                elif(hit.alignment_length < longest_hit_len):
                    # Only take the longest hits.
                    continue
                elif (curr_designation > max_designation):
                    # Do not take poorly designated seqs
                    continue
                elif(lowest_designation > curr_designation):
                    # Throw away any previous hits which have a designation
                    # larger than curr_designation.
                    gen_locs[:] = []

                # We accept this hit (for now).
                longest_hit_len = hit.alignment_length
                hit.designation = hit.calculate_designation_full(curr_designation)
                gen_locs.append(hit)
                lowest_designation = curr_designation

            # Remove the query sequence from the hit to save memory.
            hit.remove_query_sequence()

        # At this point, the designation of all members of gen_locs are either:
        # a) designated as only 'threes', only 'fours', only 'fives', and so on, or
        # b) either 'ones' or 'twos' but not yet designated as either.
        # The next section actually records those designations.
        #
        # We also need to handle seqs without any worthy hits.

        if len(gen_locs) == 0:
            # No genomic locations, either due to prost-post-filtering or due to
            # truly not hitting to the genome.  Store any of the
            # prost-post-filtered hits to each short_seqs "no_hit_gen_locs".
            # Designate as a "no_hit".
            self.designate_sequence(query_sequence, DESIGNATION_NO_HIT,
                    max_locations_allowed, None, no_hit_gen_locs)
        elif gen_locs[0].designation is None:
            # Designation: 1 or 2 (for all members of gen_locs[]).
            self.designate_sequence(query_sequence, DESIGNATION_ONE_OR_TWO,
                    max_locations_allowed, gen_locs)
        else:
            # Designation: 3 or greater (for all members of gen_locs[]).
            self.designate_sequence(query_sequence, lowest_designation,
                    max_locations_allowed, gen_locs)

    def normalize(self):
        """Normalize per-sample counts for each ShortSeq.

        Note:
            This method performs the following duties (and caches these
            results in various places):

                1.  Calculates per-sample totals at the ShortSeqs level.
                2.  Calculates per-sample normalization of each ShortSeq
                    object's per-sample read totals w/r/t the ShortSeqs
                    singleton's per-sample totals.

        """
        # Note: self (ie ShortSeqs) cannot be empty here, so we do not check if
        # it is empty.  It cannot be empty because Prost bails if search.fa is
        # empty earlier on.
        repr_short_seq = next(self.itervalues())
        # get number of samples
        num_samples = repr_short_seq.samples_counts.num_samples

        # Use a SamplesCounts object to store the total for each sample
        self.per_sample_read_totals = SamplesCounts(num_samples)

        # Calc & cache per-sample read count totals across ShortSeqs
        # singleton (ignoring sequences that do not hit to the genome).
        progress = Progress("Normalization: calculating per-sample totals.")
        for short_seq in self.itervalues():
            if short_seq.is_no_hit:
                # Skip no hit seqs
                continue
            self.per_sample_read_totals += short_seq.samples_counts
        progress.done()

        # Normalize each ShortSeq's per-sample calcs.
        progress = Progress("Normalization: normalizing read counts.")
        for short_seq in self.itervalues():
            if short_seq.is_no_hit:
                # Skip no hit seqs
                continue
            short_seq.samples_counts.normalize_samples(self.per_sample_read_totals)
        progress.done()

    def designation_step_two(self):
        """Designates each ShortSeq object as either a 'one' or 'two'.

        This stage of designation only examines ShortSeq objects that had been
        previously designated DESIGNATION_ONE_OR_TWO, and then designates each
        of those as either DESIGNATION_ONE or DESIGNATION_TWO.

        Internally, this method iterates over the list of ShortSeq objects
        sorted by maximum normalized count, so that if two perfect hit sequences
        are within_wiggle of one another, the ShortSeq with the larger sum_norm
        is designated as the 'one', which the other is designated as a 'two'.

        """

        # All the ones and twos
        ones_and_twos = (short_seq for short_seq in self.itervalues() if
                short_seq.designation_integer == DESIGNATION_ONE_OR_TWO)

        # All the ones and twos sorted by max normalized count.
        ones_and_twos_sorted_by_sum_norm = sorted(ones_and_twos,
                key=operator.attrgetter('sum_norm', 'seq_str'), reverse=True)

        progress = Progress(
                "Designation step TWO",
                100, len(ones_and_twos_sorted_by_sum_norm))
        for short_seq in ones_and_twos_sorted_by_sum_norm:
            progress.progress()
            if (not self.first_sequence_designated_one_has_been_designated()):
                # Skip wiggle calculation.
                # Special case for the first sequence designated.
                self.designate_sequence(short_seq.seq_str, DESIGNATION_ONE)
            else:
                # otherwise check if it's a two or a one
                if (self.is_within_wiggle_of_ones(short_seq)):
                    # Yes, there is at least one 'one' sequence within_wiggle of
                    # gen_locs.  This makes this a 'two'.
                    self.designate_sequence(short_seq.seq_str, DESIGNATION_TWO)
                else:
                    # Designation: 1 (it wasn't a two).
                    self.designate_sequence(short_seq.seq_str, DESIGNATION_ONE)
        progress.done()

    def load_from_cache(self):
        ### NOTE - Currently broken! We stopped using GenomicLocation and instead are now just using AlignmentExecutionHits.
        """Loads genomic locations and sequence designations from the database
        cache.

        Assumes that ShortSeqs has already been populated with minimal ShortSeq
        objects from the sample input files.
        """

        # 1. Load up each sequence from ShortSeqs.
        # 2. For exah sequence, try to find it in the db cache.
        # 3. If successful, load up each genomic_location from the db cache.
        # 4. Build up gen_locs, then designate the sequence.
        progress = Progress(
                "Reading from database cache, this may take some time",
                100, len(self))
        for short_seq in self.itervalues():
            progress.progress()
            try:
                seq_cached = ShortSeqCache.get(
                    ShortSeqCache.sequence == short_seq.seq_str)
            except DoesNotExist:
                continue
            seq_str = seq_cached.sequence
            short_seq = self[seq_cached.sequence]
            designation_integer = seq_cached.designation_integer
            gen_locs = []

            gen_locs_in_cache = (GenomicLocationCache.
                    select().
                    join(ShortSeqCache).
                    where(ShortSeqCache.sequence == seq_str))

            for gl in gen_locs_in_cache:
                if gl.designation:
                    ### NOTE -  Currently broken! We stopped using GenomicLocation and instead are now just using AlignmentExecutionHits.
                    gen_locs.append(GenomicLocationWithDesignation(
                        gl.lg, gl.start, gl.end, gl.designation))
                else:
                    ### NOTE -  Currently broken! We stopped using GenomicLocation and instead are now just using AlignmentExecutionHits.
                    gen_locs.append(GenomicLocation(
                        gl.lg, gl.start, gl.end))
            self.designate_sequence(seq_str, designation_integer, gen_locs)
            short_seq.cached = True
        progress.done()


    def perform_annotations(self, alignments, species):
        """... perform annotations ...

        Arguments:
            alignments (Alignments): The singleton Alignments object.
            species (str): The species under investigation.
        """

        for annotation_alignment in alignments.annotation_alignments:
            annotation_cls = annotation_alignment.annotation_cls
            if annotation_cls == MirbaseMirReverseAnnotation:
                self._perform_reverse_annotation(annotation_alignment, species)
            else:
                self._perform_annotation(annotation_alignment, species)

    def _perform_annotation(self, annotation_alignment, species):
        """... perform a single annotations (eg only a single BBMap search) ...

        NOTE: Currently, annotation hits are required to be full_length.  We
            could change this to be configurable in the future if desired.

        Arguments:
            annotation_alignment: A single AnnotationAlignment.
            species (str): The species under investigation.
        """

        annotation_cls = annotation_alignment.annotation_cls
        for hit in annotation_alignment.alignment_execution.hits_filtered():
            # Note - Here we have *hard-coded* the fact that an alignment to an
            # annotation sequence needs to be full length.
            if (    (not hit.is_no_hit) and
                    hit.is_full_length and
                    (hit.reference_start < hit.reference_end)):
                try:
                    short_seq = self[hit.query_sequence]
                except KeyError:
                    # In cases where we are using alignment results that were
                    # generated outside of the current execution of Prost
                    # (perhaps with less stringent min seq count), we may find
                    # hits whose query sequences have been removed from
                    # ShortSeqs due to low seq count.  In that case, simply
                    # continue.
                    continue

                if issubclass(annotation_cls, MirbaseAnnotation):
                    annotation = annotation_cls(hit, species)
                    if annotation not in short_seq.annotations:
                        short_seq.annotations.append(annotation)
                elif issubclass(annotation_cls, BiomartOtherRNAAnnotation):
                    # Note that, like the other annotation hits, we create
                    # a single BiomartOtherRNAAnnotation per hit; however, the
                    # uncompressed output is spread across two columns (first
                    # column is other_ncRNA annotation (e.g. Mir99b), while
                    # second column is the biotype (of the first column) (e.g.
                    # lincRNA, miRNA, snoRNA, etc).
                    annotation = annotation_cls(hit)
                    if annotation not in short_seq.annotations:
                        short_seq.annotations.append(annotation)

            # Remove the query sequence from the hit to save memory.
            hit.remove_query_sequence()

    def _perform_reverse_annotation(self, annotation_alignment, species):
        """... performs the REVERSE annotation ...

        This adds an extra column representing the reverse miR annotation, in
        which the query sequence is miRBase mature miRs, and the reference
        sequences is our compressed reads.

        If a normal "forward" annotation is present, then *no* reverse
        annotation will be performed.

        Arguments:
            annotation_alignment: A single AnnotationAlignment.
            species (str): The species under investigation.
        """


        annotation_cls = annotation_alignment.annotation_cls
        assert(annotation_cls == MirbaseMirReverseAnnotation)
        for hit in annotation_alignment.alignment_execution.hits_full_len_100_perc_core():

            # Note - Above we have *hard-coded* the fact that a reverse alignment
            # a full length perfect match.
            if (hit.reference_start < hit.reference_end):

                try:
                    short_seq = self[hit.reference_sequence_name]
                except KeyError:
                    # In cases where we are using alignment results that were
                    # generated outside of the current execution of Prost
                    # (perhaps with less stringent min seq count), we may find
                    # hits whose query sequences have been removed from
                    # ShortSeqs due to low seq count.  In that case, simply
                    # continue.
                    continue

                if any(a.__class__ == MirbaseMirAnnotation and
                        a.kind == a.Kind.species for a in short_seq.annotations):
                    # This short_seq already has a forward annotation, do not
                    # give it a reverse annotation.
                    # Preventing rev_annos when fwd_annos exist in gen_loc_bins
                    # is implemented elsewhere, in mirbase_mir_reverse_annotations().
                    
                    continue

                annotation = annotation_cls(hit, species)
                if annotation not in short_seq.annotations:
                    short_seq.annotations.append(annotation)

            # Remove the query sequence from the hit to save memory.
            hit.remove_query_sequence()


class Bins(object):
    """Abstract Base Class for Bins."""
    __metaclass__ = ABCMeta

    def __len__(self):
        return len(self._bins)

    def __init__(self, short_seqs_singleton):
        """
        Args:
            short_seqs_singleton (ShortSeqs): The ShortSeqs singleton.

        """
        self.samples_counts_totals = None
        self._short_seqs = short_seqs_singleton
        self._bins = set()
        self._bins_indexed_by_seq = {}
        self._debug_seq_binned_count = 0
        self._current_idx = 1
        self.per_sample_read_totals = None

    def __iter__(self):
        """Make Bins iterable by returning the _bins iterator."""
        return self._bins.__iter__()

    def next(self):
        """Make Bins iterable by returing result of _bins next() method."""
        return self._bins.next()

    @abstractmethod
    def perform_binning(self):
        """Perform the Binning."""
        pass

    def manipulate_per_sample_counts(self):
        """Perform various manipulations on the per-sample counts found in this
        Bins object, its child Bin objects, and those child Bin objects' child
        ShortSeq objects.

        Note:
            This method performs the following duties (and caches these
            results in various places):

                1.  Calculates per-sample totals at the Bin and Bins levels.
                2.  Calculates per-sample normalization of each Bin's per-sample
                    read totals w/r/t this Bins object's per-sample totals.

        """
        self._calculate_per_sample_totals()
        self._normalize()

    def _calculate_per_sample_totals(self):

        # Only proceed if there are bins to work on.
        if self._bins:

            # get number of samples
            repr_bn = next(iter(self._bins))
            repr_seq = self._short_seqs[repr_bn.main_short_seq_str]
            num_samples = repr_seq.samples_counts.num_samples

            # Use a SamplesCounts object to store the total for each sample
            self.per_sample_read_totals = SamplesCounts(num_samples)

            # Tally-ho
            for bn in self._bins:
                bn.manipulate_per_sample_counts(self._short_seqs)
                self.per_sample_read_totals += bn.per_sample_read_totals

    def _normalize(self):
        for bn in self._bins:
            # normalize bin read counts per-sample
            bn.per_sample_read_totals.normalize_samples(self.per_sample_read_totals)

    def _start_bin(self, short_seq, bn):
        """Creates (starts) a new bin with a ShortSeq obj as its first member.

        Note:
            Actually, Bins only does accounting on bn.  The subclass is
            technically responsible for creating a new Bin object, and then
            delegates the housekeeping to Bins.

        Args:
            short_seq (ShortSeq): The ShortSeq starting the bin.
            bn (Bin): The subclass of Bin object created by the subclass of Bins.

        Returns:
            Bin: The bin that was started.

        """
        assert self.get_bin_by_seq(short_seq.seq_str) == None
        self._debug_seq_binned_count += 1
        self._current_idx += 1
        self._bins.add(bn)
        self._bins_indexed_by_seq[short_seq.seq_str] = bn
        return bn

    def _add_to_bin(self, main_short_seq_str, joining_short_seq_str):
        """Adds a ShortSeq to an existing bin (with seq_strs).

        Args:
            main_short_seq_str (str): The main ShortSeq of this bin. Different
                subclasses of bins affix different meaning to what 'main' means.
                For example, for a GenLocBin, the main short_seq is the
                short_seq that started the bin.  See relevant docs in each of
                the subclasses of Bins.
            joining_short_seq_str (str): The ShortSeq to add to the bin.

        """
        assert self.get_bin_by_seq(joining_short_seq_str) == None
        self._debug_seq_binned_count += 1
        bn = self.get_bin_by_seq(main_short_seq_str)
        bn.append(joining_short_seq_str)
        self._bins_indexed_by_seq[joining_short_seq_str] = bn

    def all_bins_sorted_by_starting_seq(self):
        return sorted(self._bins, key=operator.attrgetter('main_short_seq_str'))

    def get_bin_by_seq(self, seq_str):
        """Retrieve (if any) the bin to which this short_seq belongs.

        Args:
            seq_str (str): The sequence string of a short_seq.

        Returns:
            Bin: The Bin containing the ShortSeq with seq_str.

        """
        return self._bins_indexed_by_seq.get(seq_str, None)

    def reset_bins_main_short_seq_to_sum_norm_prefer_low_designation(self):
        """Walk through each bin and set its main_short_seq to its member with
        the highest sum_norm, preferring lower designations over others.

        Used currently by SeedBins and AnnotationsBins.

        """
        for bn in self._bins:
            bn.reset_bin_main_short_seq_to_sum_norm_prefer_low_designation(self._short_seqs)


class SeedBins(Bins):
    """Container for seed region bins."""

    _idx_prefix = 's'

    def __init__(self, short_seqs_singleton):
        """
        Args:
            short_seqs_singleton (ShortSeqs): The ShortSeqs singleton.

        """
        super(SeedBins, self).__init__(short_seqs_singleton)
        self._bins_indexed_by_seed = {}

    def __repr__(self):
        asdf = {'idx': pprint.pformat(self._bins_indexed_by_seed),
                'lst': pprint.pformat(self._bins)}
        return pprint.pformat(asdf)

    def perform_binning(self, annotation_bins):
        """Perform the Binning by seed region."""

        for short_seq in self._short_seqs.sorted_by_designation_then_sum_norm():
            if short_seq.is_no_hit:
                continue
            anno_bin = annotation_bins.get_bin_by_seq(short_seq.seq_str)
            seed_bn = self._bins_indexed_by_seed.get(short_seq.seed, None)
            if anno_bin is None:
                # We do not include short seqs without annotation in seed bins
                pass
            elif seed_bn:
                main_short_seq_str = seed_bn.main_short_seq_str
                self._add_to_bin(main_short_seq_str, short_seq.seq_str)
            else:
                self._start_bin(short_seq)

        # Choose the highest sum_norm seq from each bin to be each bin's
        # main_short_seq preferring low designations over others.
        self.reset_bins_main_short_seq_to_sum_norm_prefer_low_designation()

    def _start_bin(self, short_seq):
        """Creates (starts) a new bin with a ShortSeq obj as its first member.

        Args:
            short_seq (ShortSeq): The ShortSeq starting the bin.

        Returns:
            Bin: The bin that was started.

        """
        bn = SeedBin(short_seq.seq_str, self._idx_prefix + str(self._current_idx))
        super(SeedBins, self)._start_bin(short_seq, bn)
        seed = short_seq.seed
        self._bins_indexed_by_seed[seed] = bn
        return bn


class AnnotationBins(Bins):
    """Container for annotation bins.

    This Bins subclass is a little bit different from the other Bins subclasses
    (e.g. SeedBins) in that it takes as input another Bins subclass
    singleton object, namely GenLocBins, and uses that as a starting point for
    annotation.

    This is somewhat of a shortcut, but also is fundamentally different than
    walking through each ShortSeq independently (see below).

    Specifically, AnnotationBins effectively re-bins the GenLocBins by
    collecting and concatenating GenLocBins which share the same annotation.
    Note that this gives a slightly different result than if we were to
    perform_binning on the raw set of ShortSeqs.  The main (?) case which is
    different is when a #1 does not have an annotation, but its bin members
    do."""

    _idx_prefix = 'a'

    def __init__(self, short_seqs_singleton, gen_loc_bins):
        """
        Args:
            short_seqs_singleton (ShortSeqs): The ShortSeqs singleton.
            gen_loc_bins (GenLocBins): The GenLocBins singleton.

        """
        super(AnnotationBins, self).__init__(short_seqs_singleton)
        self._gen_loc_bins = gen_loc_bins
        self._bins_indexed_by_annotations_names = {}
        self.arm_5p_3p_pairs = []
        self._debug_seq_attempted_binned_count = 0

    def _start_bin(self, gen_loc_bin, annotations):
        """Creates (starts) a new bin with a ShortSeq obj as its first member.

        Args:
            gen_loc_bin (GenLocBin): The GenLocBin starting the bin.
            annotations ((Annotation,)): The sorted Annotation tuple starting this bin.

        Note:
            Currently we are not allowing *any* no_hit_seqs inside of bins, and
            hence we don't even need to provide logic to deal with no_hit_seqs
            in here (since they can't even get into gen_loc_bins, then since
            annotation_bins are created from gen_loc_bins, it's impossible for
            them to get in here).  If in the future we need to allow no_hit_seqs
            in some/all bins, then messy/smelly logic will be required here.

        Returns:
            Bin: The bin that was started.

        """
        main_short_seq_str = gen_loc_bin.main_short_seq_str
        main_short_seq = self._short_seqs[main_short_seq_str]
        bn = AnnotationBin(main_short_seq_str,
                    self._idx_prefix + str(self._current_idx),
                    annotations)
        super(AnnotationBins, self)._start_bin(main_short_seq, bn)
        annotations_names = tuple(a.name for a in annotations)
        self._bins_indexed_by_annotations_names[annotations_names] = bn

        # And of course, add gen_loc_bin's joining members as well
        for short_seq_str in gen_loc_bin.joining_short_seq_strs:
            super(AnnotationBins, self)._add_to_bin(
                main_short_seq_str, short_seq_str)
        return bn

    def _add_to_bin(self, gen_loc_bin, annotation_bin):
        """"Add" all the ShortSeq members in gen_loc_bin to annotation_bin.

            Args:
                gen_loc_bin (GenLocBin): The GenLocBin we are appending to bn.
                annotation_bin (AnnotationBin): The AnnotationBin to which we
                    are appending the members of GenLocBin.
        """
        main_short_seq_str = annotation_bin.main_short_seq_str
        for short_seq_str in gen_loc_bin.members:
            super(AnnotationBins, self)._add_to_bin(
                main_short_seq_str, short_seq_str)

    def perform_binning(self):
        """Perform the Binning by annotation.

        Note:
            Currently we are not allowing *any* no_hit_seqs inside of *any* kind
            of bins.  Nonetheless, providing a failsafe inside just in case.
            See _start_bin().
        """
        # Sort by short_seq_str for reproducibility.

        # First pass - start with forward annotations only.
        for gen_loc_bin in sorted(self._gen_loc_bins,
                        key=operator.attrgetter('main_short_seq_str')):

            if self._short_seqs[gen_loc_bin.main_short_seq_str].is_no_hit:
                # Just in case, adding a short circuit here just in case.
                continue

            species_mir_annos = \
                gen_loc_bin.mirbase_mir_annotations(self._short_seqs)[SPECIES_IDX]
            species_mir_annos_names = tuple(a.name for a in species_mir_annos)

            if species_mir_annos == ():
                # Doesn't have normal in_species forward annotations. Skip.
                continue

            bn = self._bins_indexed_by_annotations_names.get(
                    species_mir_annos_names, None)
            if bn:
                self._add_to_bin(gen_loc_bin, bn)
            else:
                self._start_bin(gen_loc_bin, species_mir_annos)

        # Second pass - reverse annotations time
        for gen_loc_bin in sorted(self._gen_loc_bins,
                        key=operator.attrgetter('main_short_seq_str')):

            if self._short_seqs[gen_loc_bin.main_short_seq_str].is_no_hit:
                # Just in case, adding a short circuit here just in case.
                continue

            species_mir_annos = \
                gen_loc_bin.mirbase_mir_annotations(self._short_seqs)[SPECIES_IDX]
            species_mir_rev_annos = \
                gen_loc_bin.mirbase_mir_reverse_annotations(self._short_seqs)[SPECIES_IDX]
            species_mir_rev_annos_names = tuple(a.name for a in species_mir_rev_annos)

            if species_mir_annos != ():
                # Has normal in_species forward annotations. Already done in
                # first pass.  Skip.
                continue
            elif species_mir_rev_annos == ():
                # Has neither forward or reverse in_species annos.  Skip.
                continue
            else:
                # Doesn't have foward, but does have reverse in_species annos. Use that.
                pass

            bn = self._bins_indexed_by_annotations_names.get(
                    species_mir_rev_annos_names, None)
            if bn:
                self._add_to_bin(gen_loc_bin, bn)
            else:
                print("Of note: Found a rev_anno_only bin_starter! gen_loc_bin.idx = {}, gen_loc_bin.main_short_seq_str = {}".format(gen_loc_bin.idx, gen_loc_bin.main_short_seq_str))
                self._start_bin(gen_loc_bin, species_mir_rev_annos)


        # Choose the highest sum_norm seq from each bin to be each bin's
        # main_short_seq preferring low designations over others.
        self.reset_bins_main_short_seq_to_sum_norm_prefer_low_designation()

    def post_binning_processing(self):
        """Perform various post binning processing.

        At time of writing this includes only arm switching detection.

        """
        self._arm_switch_detection()

    def _arm_switch_detection(self):
        """Perform arm switch detection.  See log_fold_changes_expr_5p_3p() for
        more detail."""

        arm_5p_annos = {}
        arm_3p_annos = {}
        for (annos_names, bn) in self._bins_indexed_by_annotations_names.iteritems():
            if len(annos_names) != 1:
                continue
            anno_name = annos_names[0]
            m = re.search(r'(.*)-([53]p$)', anno_name, re.IGNORECASE)
            if m and len(m.groups()) == 2:
                name, arm = m.groups()
                if arm == '5p':
                    arm_5p_annos[name] = bn
                elif arm == '3p':
                    arm_3p_annos[name] = bn
                else:
                    raise ControlFlowException, \
                            "ERR415: Shouldn't be possible to reach here."

        # Retrieve those that have exactly one 5p and 3p
        have_exactly_one_5p_and_3p = sorted(
            set(arm_5p_annos).intersection(arm_3p_annos))

        for species_hairpin_anno_name in have_exactly_one_5p_and_3p:
            bn_5p = arm_5p_annos[species_hairpin_anno_name]
            bn_3p = arm_3p_annos[species_hairpin_anno_name]
            arm_5p_3p_pair = (species_hairpin_anno_name, bn_5p, bn_3p)
            self.arm_5p_3p_pairs.append(arm_5p_3p_pair)


class GenLocBins(Bins):
    """Container for genomic location bins.

    Each short sequence is thrown into a single bin.  Each bin can contain
    one or more sequence strings.  Each bin is started by a "starter short
    sequence".  Bins are started roughly as follows:

    1. For sequences sorted first by best designation, then by sum_norm:
    2. If the sequence is not within_wiggle of a bin_starter, then create a bin
       for it.

    Sequences are defined as 'ambiguous' if they are within_wiggle of more than
    one bin_starter.

    Non-starter sequences are thrown into a bin if they are within_wiggle of
    exactly one bin_starter.
    """

    _idx_prefix = 'gl'

    def __init__(self, short_seqs_singleton):
        """
        Args:
            short_seqs_singleton (ShortSeqs): The ShortSeqs singleton.

        """
        super(GenLocBins, self).__init__(short_seqs_singleton)
        self._bins_indexed_by_designation = {}
        self._non_one_bin_starters = set()
        self._ambiguous = set()
        self._probation = set()
        self._debug_seq_attempted_binned_count = 0
        self.mirror_pairs = []

    def __repr__(self):
        asdf = {'idx': pprint.pformat(self._bins_indexed_by_designation),
                'lst': pprint.pformat(self._bins)}
        return pprint.pformat(asdf)

    def perform_binning(self):
        """Perform the Binning by Genomic Location"""
        if DEBUG_TIMING:
            progress = Progress('\n\t_bin_ones')
        self._bin_ones()
        if DEBUG_TIMING:
            progress.done()
            progress = Progress('\n\t_bin_others')
        self._bin_others()
        if DEBUG_TIMING:
            progress.done()
            progress = Progress('\n\t_bin_probations')
        self._bin_probations()
        if DEBUG_TIMING:
            progress.done()

    def post_binning_processing(self, conf):
        """Perform various post binning processing.

        At time of writing this includes only mirror-miRNA detection.  However,
        see the TODO below.

        TODO: Profile time taken to perform the following three
            all-bins-iterations, and determine if folding two of those
            iterations into one would speed things up at the cost of code
            complexity:

                1. calc per-bin per-sample totals for all bins
                2. read count normalization
                3. mirror-miRNA detection
        """
        # For now, set a small number for max_gen_loc_len...
        self._mirror_mir_detection(TODO__HARDCODE_SMALL_MAX_GEN_LOC_LEN_FOR_MIRROR_MIR_DETECTION)

    def _start_bin(self, short_seq):
        """Creates (starts) a new bin with a ShortSeq obj as its first member.

        Args:
            short_seq (ShortSeq): The ShortSeq starting the bin.

        Returns:
            Bin: The bin that was started.

        """
        bn = GenLocBin(short_seq.seq_str, self._idx_prefix + str(self._current_idx))
        super(GenLocBins, self)._start_bin(short_seq, bn)
        desig = short_seq.designation_integer
        self._bins_indexed_by_designation[desig] = (
                self._bins_indexed_by_designation.get(desig, []))
        self._bins_indexed_by_designation[desig].append(bn)
        if not short_seq.is_one:
            self._non_one_bin_starters.add(short_seq.seq_str)
        return bn

    def _bin_ones(self):
        """Bin the sequences designated 'ones'."""
        for short_seq in self._short_seqs.ones():
            self._debug_seq_attempted_binned_count += 1
            self._start_bin(short_seq)

    def _bin_others(self):
        """

        Algorithm is roughly as follows:
            Walking in order of best designation then sum_norm:
            1. mark ambiguous if within_wiggle of more than one 'one'.
            2. if within_wiggle of just one 'one', put it in that bin
            3. if within_wiggle of no ones, *might* be a bin starter
            4. mark ambiguous if within_wiggle of more than one bin_starter
            5. if within_wiggle of just one bin_starter, mark it probation
            6. if within_wiggle of no bin_starters, it's a bin starter
        """
        # Dec/13/13 - Thomas confirms that this is the right sort order.
        for short_seq in \
                self._short_seqs.non_ones_sorted_by_designation_then_sum_norm():

            if short_seq.is_no_hit:
                # Stat: DESIGNATION_NO_HIT, do not bin
                continue

            if short_seq.max_locations_allowed_exceeded:
                # Stat: Too many locations to store, do not bin
                continue

            self._debug_seq_attempted_binned_count += 1

            # lists of seqs for which this seq is within_wiggle
            wigglers = self._short_seqs.within_wiggle_of(short_seq, True)
            one_wigglers = [s.seq_str for s in wigglers if s.is_one]

            if (len(one_wigglers) == 1):
                # Stat: bin with the ones
                one = self._short_seqs[one_wigglers[0]]
                self._add_to_bin(one.seq_str, short_seq.seq_str)
            elif(len(one_wigglers) > 1):
                # Stat: ambiguous
                # TODO: next two lines should be a single line (single point of
                # responsibility)
                short_seq.set_ambiguous(one_wigglers)
                self._ambiguous.add(short_seq.seq_str)
            else:
                # Look for bin_starters (i.e. non_ones) which are within_wiggle
                # of short_seq.

                wigglers = self._short_seqs.within_wiggle_of(short_seq, True)
                non_one_wigglers = {s.seq_str for s in wigglers if not s.is_one}
                non_one_bin_starters_within_wiggle = \
                        non_one_wigglers.intersection(self._non_one_bin_starters)
                if (len(non_one_bin_starters_within_wiggle) > 0):
                    # Stat: probation
                    self._probation.add(short_seq.seq_str)
                else:
                    # Stat: non_one_bin_starter
                    self._start_bin(short_seq)

    def _bin_probations(self):
        """Bin the sequences previously on probation (if possible)."""

        for short_seq in (self._short_seqs[s] for s in self._probation):
            wigglers = self._short_seqs.within_wiggle_of(short_seq, True)
            non_one_wigglers = {s.seq_str for s in wigglers if not s.is_one}
            non_one_bin_starters_within_wiggle = list(
                    non_one_wigglers.intersection(self._non_one_bin_starters))
            if(len(non_one_bin_starters_within_wiggle) > 1):
                # Stat: probation becomes ambiguous
                # TODO: next two lines should be a single line (single point of
                # responsibility)
                short_seq.set_ambiguous(non_one_bin_starters_within_wiggle)
                self._ambiguous.add(short_seq.seq_str)
                self._debug_seq_binned_count += 1
            else:
                # Stat: non_one_bin_starter
                bin_starter = \
                        self._short_seqs[non_one_bin_starters_within_wiggle[0]]
                self._add_to_bin(bin_starter.seq_str, short_seq.seq_str)

    def _calculate_per_sample_totals(self):
        """Calculate per-sample read totals."""

        # Only proceed if there are bins to work on.
        if self._bins:

            # get number of samples
            repr_bn = next(iter(self._bins))
            repr_seq = self._short_seqs[repr_bn.main_short_seq_str]
            num_samples = repr_seq.samples_counts.num_samples

            # Use a SamplesCounts object to store the total for each sample
            self.per_sample_read_totals = SamplesCounts(num_samples)

            # Tally-ho
            for bn in self._bins:
                bn.manipulate_per_sample_counts(self._short_seqs)
                main_seq = self._short_seqs[bn.main_short_seq_str]
                if main_seq.is_no_hit:
                    # Skip no hit seqs
                    continue
                self.per_sample_read_totals += bn.per_sample_read_totals

    def _mirror_mir_detection(self, max_gen_loc_len):
        """Mirror miR detection module.

        A first pass at mirror miR detection.  Better mirror-miRNA detection
        will be available after novel miR detection is implemented.

        Args:
            max_gen_loc_len (int): The maximum number of genomic locations that
                a bin to be included in the detection module.
        """
        mirror_bins = {}

        # Create the bins
        for bn in self._bins:
            starter_seq = self._short_seqs[bn.main_short_seq_str]
            if (len(starter_seq.genomic_locations) > max_gen_loc_len):
                continue
            for gen_loc in starter_seq.genomic_locations:
                mirror_fp = to_mirror_fingerprint(gen_loc)
                mirror_bins[mirror_fp] = mirror_bins.get(mirror_fp, [])
                mirror_bins[mirror_fp].append((bn, gen_loc))


        # Find the potential mirror miRs
        for mirror_bin in mirror_bins.itervalues():
            for i in xrange(len(mirror_bin)):
                for j in xrange(i+1, len(mirror_bin)):
                    (bn1, gen_loc1) = mirror_bin[i]
                    (bn2, gen_loc2) = mirror_bin[j]
                    if bn1 == bn2:
                        continue
                    if gen_loc1.within_mirror(gen_loc2):
                        # Always put the pair with the "lowest" gen_loc first
                        # for consistency
                        if gen_loc1 < gen_loc2:
                            mirror_pair = (bn1, gen_loc1, bn2, gen_loc2)
                        else:
                            mirror_pair = (bn2, gen_loc2, bn1, gen_loc1)
                        self.mirror_pairs.append(mirror_pair)

        # sort for consistency
        self.mirror_pairs = sorted(self.mirror_pairs, key=
                                lambda(mp): (mp[1], mp[3], mp[0].idx, mp[2].idx))


class Bin(SlotPickleMixin):
    """A generic Bin object, used to group together short_seq objects according
    to various criteria.  At time of writing, this includes binning by genomic
    location, by seed, and by annotation.

    Each short sequence is thrown into a single bin.  Each bin can contain
    one or more sequence strings.  Each bin is started by a "starter short
    sequence", which may or may not have meaning (for example, the bin starter
    in the genomic location bins is generally the most expressed region within a
    small wiggle region, whereas the bin starter in the seed bins will likely be
    just the first short_seq encountered with that bin's seed).

    A 'Bin' object contains a starter sequence, and zero or more joining
    members.

    The class is not iterable itself to ensure that iteration over members
    remains explicit (that is, are you iterating over all members, or just
    joining members?).

    Properties:
        members: The sequence strings of ShortSeq objects that are members of
            this bin.
    """
    __slots__ = (
            'main_short_seq_str',
            'joining_short_seq_strs',
            'idx',
            'per_sample_read_totals',

            # Alternatively cut on 5prime (templated or untemplated?)
            'per_sample_seed_shifted_totals',
            # NTs 2-8 -- ins/del/mismatch
            'per_sample_seed_edited_totals',
            # NTs 9-12
            'per_sample_central_edited_totals',
            # NTs 13-17
            'per_sample_3p_supplementary_edited_totals',

            # Shorter/Longer/Untemplated-addition

            # Alternatively cut on 3p end
            'per_sample_3p_alt_cut_totals',
            # 3p mismatch, i.e. putative untemplated addition
            'per_sample_3p_mismatch_totals',
            # NTs 1, 9-12 and/or 17 onwards
            'per_sample_other_edited_totals',
            )

    @property
    def members(self):
        """The members of this bin (as sequence strings)."""
        return [self.main_short_seq_str] + self.joining_short_seq_strs

    def __init__(self, main_short_seq_str, idx):
        """Initialization.

        Args:
            main_short_seq_str (str): The sequence string of the main ShortSeq
                of this Bin.
            joining_short_seq_strs ([str]): List of the sequence strings of 'ShortSeq's
                that joined this bin.
        """
        self.main_short_seq_str = main_short_seq_str
        self.joining_short_seq_strs = []
        self.idx = idx
        self.per_sample_read_totals = None
        self.per_sample_seed_shifted_totals = None
        self.per_sample_seed_edited_totals = None
        self.per_sample_central_edited_totals = None
        self.per_sample_3p_supplementary_edited_totals = None
        self.per_sample_3p_alt_cut_totals = None
        self.per_sample_3p_mismatch_totals = None
        self.per_sample_other_edited_totals = None

    def has_indelnts(self, short_seqs):
        """Do any sequences in this Bin have indels?

        Arguments:
            short_seqs: The ShortSeqs singleton.
        """
        return any([short_seqs[ss].has_indelnts for ss in self.members])

    def append(self, joining_short_seq_str):
        """Append a joining member ShortSeq to this Bin."""
        self.joining_short_seq_strs.append(joining_short_seq_str)

    def manipulate_per_sample_counts(self, short_seqs):
        """Perform various manipulations on the per-sample counts found in this
        Bin object and its member ShortSeq objects.

        Note:
            This method performs the following duties (and caches these results
            in various places):

                1.  Calculates per-sample totals across this Bin.

        Args:
            short_seqs (ShortSeqs): The ShortSeqs singleton.

        """
        main_seq = short_seqs[self.main_short_seq_str]
        num_samples = main_seq.samples_counts.num_samples
        self.per_sample_read_totals = SamplesCountsWithNorms(num_samples)

        for short_seq in (short_seqs[m] for m in self.members):

            if short_seq.is_no_hit:
                # Skip no hit seqs
                continue

            # Calc & cache per-sample read count totals across this Bin
            self.per_sample_read_totals += short_seq.samples_counts

    def per_sample_perc_seed_shifted(self, na):
        """Calculate the per-sample percentage of this Bin's members which have
        their seed shifted with respect to the Bin's main sequence.

        Args:
            na (bool): If na is true, then instead return a list of same length
                but filled with 'NA's. Useful for output purposes.

        Returns:
            [float]: An array of floats representing the per-sample percentage
                of this Bin's members which have been seed shifted.  Or an array
                of 'NA's if 'na' is True.

        """
        if na:
            return self.per_sample_seed_shifted_totals.na()
        else:
            return self.per_sample_seed_shifted_totals.perc(self.per_sample_read_totals)

    def per_sample_perc_seed_edited(self, na):
        """Calculate the per-sample percentage of this Bin's members which have
        their seed edited with respect to the Bin's main sequence.

        Args:
            na (bool): If na is true, then instead return a list of same length
                but filled with 'NA's. Useful for output purposes.

        Returns:
            [float]: An array of floats representing the per-sample percentage
                of this Bin's members which have been seed edited.  Or an array
                of 'NA's if 'na' is True.

        """
        if na:
            return self.per_sample_seed_edited_totals.na()
        else:
            return self.per_sample_seed_edited_totals.perc(self.per_sample_read_totals)

    def per_sample_perc_central_edited(self, na):
        """Calculate the per-sample percentage of this Bin's members which have
        their central region edited with respect to the Bin's main sequence.

        Args:
            na (bool): If na is true, then instead return a list of same length
                but filled with 'NA's. Useful for output purposes.

        Returns:
            [float]: An array of floats representing the per-sample percentage
                of this Bin's members which have been central region
                edited.  Or an array of 'NA's if 'na' is True.

        """
        if na:
            return self.per_sample_central_edited_totals.na()
        else:
            return self.per_sample_central_edited_totals.perc(
                self.per_sample_read_totals)

    def per_sample_perc_3p_supplementary_edited(self, na):
        """Calculate the per-sample percentage of this Bin's members which have
        their supplementary region edited with respect to the Bin's main sequence.

        Args:
            na (bool): If na is true, then instead return a list of same length
                but filled with 'NA's. Useful for output purposes.

        Returns:
            [float]: An array of floats representing the per-sample percentage
                of this Bin's members which have been supplementary region
                edited.  Or an array of 'NA's if 'na' is True.

        """
        if na:
            return self.per_sample_3p_supplementary_edited_totals.na()
        else:
            return self.per_sample_3p_supplementary_edited_totals.perc(
                self.per_sample_read_totals)

    def per_sample_perc_3p_alt_cut(self, na):
        """Calculate the per-sample percentage of this Bin's members which have
        been alternatively cut on their 3p end with respect to the Bin's main
        sequence.

        Args:
            na (bool): If na is true, then instead return a list of same length
                but filled with 'NA's. Useful for output purposes.

        Returns:
            [float]: An array of floats representing the per-sample percentage
                of this Bin's members which have been 3p alternatively cut. Or
                an array of 'NA's if 'na' is True.

        """
        if na:
            return self.per_sample_3p_alt_cut_totals.na()
        else:
            return self.per_sample_3p_alt_cut_totals.perc(
                self.per_sample_read_totals)

    def per_sample_perc_3p_mismatch(self, na):
        """Calculate the per-sample percentage of this Bin's members which have
        3p mismatches (putative untemplated addition) with respect to the Bin's
        main sequence.

        Args:
            na (bool): If na is true, then instead return a list of same length
                but filled with 'NA's. Useful for output purposes.

        Returns:
            [float]: An array of floats representing the per-sample percentage
                of this Bin's members which have 3p mismatches (putatitive
                untemplated addition). Or an array of 'NA's if 'na' is True.

        """
        if na:
            return self.per_sample_3p_mismatch_totals.na()
        else:
            return self.per_sample_3p_mismatch_totals.perc(
                self.per_sample_read_totals)

    def per_sample_perc_other_edited(self, na):
        """Calculate the per-sample percentage of this Bin's members which have
        nts 1, 9-12, or 17 onwards edited with respect to the Bin's main sequence.

        Args:
            na (bool): If na is true, then instead return a list of same length
                but filled with 'NA's. Useful for output purposes.

        Returns:
            [float]: An array of floats representing the per-sample percentage
                of this Bin's members which have been "other edited" edited.  Or
                an array of 'NA's if 'na' is True.

        """
        if na:
            return self.per_sample_other_edited_totals.na()
        else:
            return self.per_sample_other_edited_totals.perc(
                self.per_sample_read_totals)

    def mirbase_mir_annotations(self, short_seqs):
        """See self._mirbase_annotations()."""
        return self._mirbase_annotations(short_seqs, MirbaseMirAnnotation)

    def mirbase_mir_reverse_annotations(self, short_seqs):
        """See self._mirbase_annotations().

        Note:
            If there exists an in_species normal "forward" annotation, then
            there we define that there is no in_species reverse annotation.

            If there exists an other_species normal "forward" annotation, then
            there we define that there is no other_species reverse annotation.
        """
        fwd_annos = self._mirbase_annotations(short_seqs, MirbaseMirAnnotation)
        rev_annos = self._mirbase_annotations(short_seqs, MirbaseMirReverseAnnotation)
        fwd_annos_in_species, fwd_annos_other_species = fwd_annos
        rev_annos_in_species, rev_annos_other_species = rev_annos

        if len(fwd_annos_in_species) > 0:
            # there are in_species forward annotations, empty the reverse annotations
            rev_annos_in_species = ()
        if len(fwd_annos_other_species) > 0:
            # there are other_species forward annotations, empty the reverse annotations
            rev_annos_other_species = ()

        # Reconsitute the reverse annotations
        rev_annos = (rev_annos_in_species, rev_annos_other_species)

        return rev_annos

    def mirbase_hairpin_annotations(self, short_seqs):
        """See self._mirbase_annotations()."""
        return self._mirbase_annotations(short_seqs, MirbaseHairpinAnnotation)

    def _mirbase_annotations(self, short_seqs, annotation_cls):
        """
        Note:
            We follow this algorithm:
            The set of the members' annotations are collected and become
            the bin's annotations.

        Args:
            short_seqs (ShortSeqs): The ShortSeqs singleton.
            annotation_cls (type): The (sub)class of Annotation desired (e.g.
                MirbaseMirAnnotation).

        Returns:
            ((Annotation,), (Annotation,)): A 2-tuple containing a sorted tuple
            of species annotations and a sorted tuple of other_species
            annotations.
        """

        annos_list = []

        for kind in annotation_cls.Kind:
            # Retrieve main member and joining member annotations.
            all_short_seq_annos = annotations_of(
                    annotation_cls, short_seqs, kind,
                    self.main_short_seq_str, *(self.joining_short_seq_strs))

            annos = tuple(sorted(all_short_seq_annos,
                        key=operator.attrgetter('name')))
            annos_list.append(annos)

        return tuple(annos_list)

    def biomart_other_rna_annotations(self, short_seqs):
        """TODO: better intro...
        Note:
            The set of the members' annotations are collected and become
            the bin's annotations.

        Returns:
            ([str], [str]): A tuple containing a list of biomart names and a
                list of biotypes.

        """
        annos = annotations_of(BiomartOtherRNAAnnotation, short_seqs, None,
                self.main_short_seq_str, *(self.joining_short_seq_strs))
        names = []
        biotypes = []
        for anno in sorted(annos, key=operator.attrgetter('name')):
            if anno.name not in names:
                names.append(anno.name)
                biotypes.append(anno.biotype)
        return (names, biotypes)

    def reset_bin_main_short_seq_to_sum_norm_prefer_low_designation(self, short_seqs):
        """Reset this Bin's main_short_seq to its member with the highest
        sum_norm preferring lower designations over others.

        TODO: This may be inefficient (O(n^2)), profile and determine if this
            has an impact on real world data sets. See note below.

        Args:
            short_seqs (ShortSeqs): The ShortSeqs singleton.
        """
        new_main_short_seq = None
        max_sum_norm = -ONE_MILLION

        # Get the ordered set of designations in this Bin.
        desigs = set()
        for short_seq in (short_seqs[s] for s in self.members):
            desigs.add(short_seq.designation_integer)
        desigs = sorted(desigs)

        # TODO: This may be inefficient, see note above.
        for desig in desigs:
            for short_seq in (short_seqs[s] for s in self.members):
                if desig == short_seq.designation_integer:
                    if short_seq.sum_norm > max_sum_norm:
                        new_main_short_seq = short_seq
                        max_sum_norm = short_seq.sum_norm
            # Break if we found an appropriate main sequence.
            if new_main_short_seq:
                break

        # Swap in the newly found main_short_seq (if it is not already the
        # main_short_seq)
        orig_main_short_seq_str = self.main_short_seq_str
        if (new_main_short_seq.seq_str != orig_main_short_seq_str):
            new_main_short_seq_str = new_main_short_seq.seq_str
            self.joining_short_seq_strs.remove(new_main_short_seq_str)
            self.joining_short_seq_strs.append(orig_main_short_seq_str)
            self.main_short_seq_str = new_main_short_seq_str


class AnnotationBin(Bin):
    """A Bin containing ShortSeq objects grouped accordinging to their
    species_miR annotation."""

    __slots__ = ('annotations',)

    def __init__(self, main_short_seq_str, idx, annotations):
        """Initialization.

        Args:
            main_short_seq_str (str): The sequence string of this Bin's main
                short_seq.
            idx (str): The index of this bin.
            annotations ((Annotations,)): The sorted tuple of annotations that
                define this Bin.

        """
        super(AnnotationBin, self).__init__(main_short_seq_str, idx)
        self.annotations = annotations


class SeedBin(Bin):

    # Make sure we don't accidentally add new attributes.
    __slots__ = []

    def mirbase_mir_annotations(self, short_seqs, gen_loc_bins):
        """See self._mirbase_annotations()."""
        return self._mirbase_annotations(short_seqs, gen_loc_bins,
                                    MirbaseMirAnnotation)

    def mirbase_hairpin_annotations(self, short_seqs, gen_loc_bins):
        """See self._mirbase_annotations()."""
        return self._mirbase_annotations(short_seqs, gen_loc_bins,
                                    MirbaseHairpinAnnotation)

    def _mirbase_annotations(self, short_seqs, gen_loc_bins, annotation_cls):
        """Calculate the mirbase annotations for this SeedBin.

        Note:
            The algorithm used is as follows:
                Collect the annotations of each member's Genomic Location Bin's
                annotations.  The set of these annotations is the set of
                annotations of this SeedBin.
        Args:
            short_seqs (ShortSeqs): The ShortSeqs singleton.
            gen_loc_bins (GenLocBins): The GenLocBins singleton.
            annotation_cls (type): The (sub)class of Annotation desired (e.g.
                MirbaseMirAnnotation).

        Returns:
            ((Annotation,), (Annotation,)): A 2-tuple containing a
                sorted-by-name tuple of species mir annotations and a
                sorted-by-name tuple of other_species mir annotations.

        """
        annos_species = set()
        annos_other_species = set()
        for short_seq_str in self.members:
            gen_loc_bin = gen_loc_bins._bins_indexed_by_seq.get(short_seq_str,
                                                                None)
            if gen_loc_bin:
                annos = gen_loc_bin._mirbase_annotations(short_seqs,
                                            annotation_cls)
                annos_species.update(annos[0])
                annos_other_species.update(annos[1])

        annos_species = tuple(sorted(annos_species,
                key=operator.attrgetter('name')))
        annos_other_species = tuple(sorted(annos_other_species,
                key=operator.attrgetter('name')))

        return (annos_species, annos_other_species)


class GenLocBin(Bin):
    """A Bin containing ShortSeq objects grouped accordinging to their
    genomic locations."""

    # Make sure we don't accidentally add new attributes.
    __slots__ = []

    def _mirbase_annotations(self, short_seqs, annotation_cls):
        """
        Note:
            For the species_mir and species_hairpin columns, we follow
            this algorithm:
                If the DESIGNATION_ONE_OR_TWO member sequences have any
                annotations, then the set of those member sequences' annotations
                become the bin's annotations.  If they do not, then the set of
                the members' annotations (i.e. all the 3's and above) are
                collected and become the bin's annotations.

            For all other columns, we follow this algorithm:
                The set of the members' annotations are collected and become the
                bin's annotations.

        Note:
            PREVIOUSLY, the algorithm differed as follows:
                For the species_mir and species_hairpin columns, the algorithm
                previously did this:
                    If a bin's starter sequence has annotations, then those
                    annotations become the bin's annotations.  Otherwise the set
                    of the members' annotations are collected and become the
                    bin's annotations.

        Returns:
            ((Annotation,), (Annotation,)): A 2-tuple containing a
                sorted-by-name tuple of species annotations and a sorted-by-name
                tuple of other_species annotations.


        """
        annos_list = []
        desigs_ones_and_twos = [DESIGNATION_ONE, DESIGNATION_TWO]

        for kind in annotation_cls.Kind:

            # Retrieve two sets of annotations:
            #   1. the set of annotations of all bin members.
            #   2. the set of annotations of the ones and twos.
            all_short_seqs = (short_seqs[s] for s in self.members)
            ones_and_twos = [s.seq_str for s in all_short_seqs
                if s.designation_integer in desigs_ones_and_twos]
            ones_and_twos_annos = annotations_of(
                    annotation_cls, short_seqs, kind,
                    *(ones_and_twos))
            all_short_seq_annos = annotations_of(
                    annotation_cls, short_seqs, kind,
                    *(self.members))

            # Determine which annotations to use
            if (kind == annotation_cls.Kind.species):
                if ones_and_twos_annos:
                    # Ones and Twos have annotations, take them.
                    annos = ones_and_twos_annos
                else:
                    # No annotations on ones and twos. Use
                    # set of joining member annotations.
                    annos = all_short_seq_annos
            else:
                # For the other columns (OtherSpecies_miR
                # OtherSpecies_hairpin other_ncRNA
                # ncRNA_biotype) simple collect the set of
                # annotations across all bin members.
                annos = all_short_seq_annos

            annos_list.append(tuple(sorted(annos, key=operator.attrgetter('name'))))

        return tuple(annos_list)

    def mir_ambig(self, short_seqs):
        return self._ambig(short_seqs, MirbaseMirAnnotation)

    def hairpin_ambig(self, short_seqs):
        return self._ambig(short_seqs, MirbaseHairpinAnnotation)

    def _ambig(self, short_seqs, annotation_cls):
        ambig = False
        kind = annotation_cls.Kind.species

        # Retrieve starter member annotations.
        main_short_seq_annos = annotations_of(annotation_cls, short_seqs, kind,
                                        self.main_short_seq_str)

        # If there are no annotations on the main ShortSeq, it is not ambiguous.
        if not main_short_seq_annos:
            return False

        # Retrieve joining member annotations.
        joining_short_seq_annos = annotations_of(annotation_cls, short_seqs,
                                        kind, *(self.joining_short_seq_strs))

        # For this species, do the seqs in this bin have more than one unique
        # MirbaseMirAnnotation or MirbaseHairpinAnnotation?
        if (joining_short_seq_annos - main_short_seq_annos):
            ambig = True

        return ambig


    def manipulate_per_sample_counts(self, short_seqs):
        """Perform various manipulations on the per-sample counts found in this
        Bin object and its member ShortSeq objects.

        Note:
            This method performs the following duties (and caches these results
            in various places):

                1.  Calculates per-sample totals across this Bin.
                2.  Calculates per-sample totals (across this Bin) for members
                    that are alteratively cut, have additions, or have editions.

        Note:
            This method is ultimately to be used by only the
            compressed_by_genomic_location code and the
            compressed_by_genomic_location_and_annotation code.

        Arguments:
            short_seqs (ShortSeqs): The ShortSeqs singleton.

        """
        main_seq = short_seqs[self.main_short_seq_str]
        num_samples = main_seq.samples_counts.num_samples
        self.per_sample_read_totals = SamplesCountsWithNorms(num_samples)
        self.per_sample_seed_shifted_totals = SamplesCounts(num_samples)
        self.per_sample_seed_edited_totals = SamplesCounts(num_samples)
        self.per_sample_central_edited_totals = SamplesCounts(num_samples)
        self.per_sample_3p_supplementary_edited_totals = SamplesCounts(num_samples)
        self.per_sample_3p_alt_cut_totals = SamplesCounts(num_samples)
        self.per_sample_3p_mismatch_totals = SamplesCounts(num_samples)
        self.per_sample_other_edited_totals = SamplesCounts(num_samples)

        for short_seq in (short_seqs[m] for m in self.members):

            mod_things = []

            if short_seq.is_no_hit:
                # Skip no hit seqs
                continue

            # Calc & cache per-sample read count totals across this Bin
            self.per_sample_read_totals += short_seq.samples_counts

            # Calc per-sample seed_shifted, seed_edited, central edited,
            # 3' supplementary edited, 3' modification, 3' untemplated addition,
            # and other edited totals.

            # For the %-columns, we skip bins not started by a designation one.
            if (main_seq.designation_integer != DESIGNATION_ONE):
                continue

            # Sanity-check assertion.
            # main_seq and short_seq should have same number of gen_locs because
            # that's how they've ended up in the same gen_loc_bin.
            assert (len(short_seq.genomic_locations) ==
                    len(main_seq.genomic_locations)), \
                    "Different sized genomic_locations\n{}\n{}\n.".format(
                        main_seq, short_seq)

            # Identify gen_loc bins that have members with conflicting types --
            # for example, one bin may have one "seed-shifted" gen_loc and one
            # "not-seed-shifted" gen_loc.  These are problematic.  Log a warning in
            # the prost.log, and do NOT include these in the calculations.
            # See Prost issue #59 on github for more info.

            # Collect modifications for each genomic location

            for gl_main, gl_ss in itertools.izip_longest(
                    main_seq.genomic_locations, short_seq.genomic_locations):

                try:
                    mod_thing = ModificationThing(gl_main, gl_ss, main_seq, short_seq)
                    mod_things.append(mod_thing)
                except ModificationThingEncounteredNAlignment:
                    # Aligned to an 'N' (i.e. cigar has an 'M').
                    # Ignore for calculations and log.
                    pmsg(
                        "Warning: While trying to calculate % miRNA "
                        "modifications, the binstarter sequence or the member "
                        "sequence aligned to one or more 'N's in the "
                        "reference genome:\n"
                        "   main_hit:   {}\n"
                        "   member_hit: {}\n".format(
                            gl_main, gl_ss)
                    )

            # It's possible that no ModificationThings were actually created
            # (for example, the BinStarter aligns to a region with an "N").
            # In that case, we simply return without making any adjustments to
            # totals.
            if len(mod_things) == 0:
                return

            # Collect any disagreements about the properties between the
            # different genomic locations about the modification properties.
            # For example, one genomic_location bin may claim "is_3p_mismatch"
            # while another may not.

            disagreements = ModificationThing.disagreements(mod_things)

            # If there's no disagreement, we simply use the first mod_thing as
            # the representative of all mod_things.

            mod_thing = mod_things[0]

            disagreement_error_msg = """
                Warning: Found two Alignment Hits from one member
                sequence of a GenomicLocationBin with a mixture of
                {} and not {} members.  NOT including these
                counts in the %-columns calculations.
                Main seq:   {}
                Member seq: {}
                Main seq's gen_locs:   {}
                Member seq's gen_locs: {}"""

            # Note: We're now tracking any disagreements about any of the 6 modification categories.
            #       Looking at past Prost runs, the vast majority are just is_3p_alt_cut and is_3p_mismatch.
            #       The others are likely just initial development bugs that were ironed out quickly.
            #       And of course, really, logically, there should never be a disagreement about
            #       any of the "edited" categories since we're only reporting "edited" categories for
            #       gen_loc_bins with DESIGNATION_ONE binstarters.
            # $ for A in */*log ; do grep "NOT including" $A ; done | column -t | sort | uniq -c
            # 1669 is_3p_alt_cut    and  not  is_3p_alt_cut    shifted  members.  NOT  including  these
            # 2512 is_3p_mismatch   and  not  is_3p_mismatch   shifted  members.  NOT  including  these
            #    9 is_other_edited  and  not  is_other_edited  shifted  members.  NOT  including  these
            #    1 shifted          and  not  seed             shifted  members.  NOT  including  these

            ## Seed Shifted? / iso_5p?

            if 'is_seed_shifted' in disagreements:
                mod = 'is_seed_shifted'
                # Note: Logically (??) it's hard to see how there could be disagreements
                #       between different genomic locations w/r/t is_seed_shifted.
                #       Nonetheless, logging is_seed_shifted disagreements just like the
                #       others.  See counts listed in note above.
                pmsg(disagreement_error_msg.format(mod, mod, main_seq,
                        short_seq, main_seq.genomic_locations,
                        short_seq.genomic_locations))
            elif mod_thing.is_seed_shifted:
                self.per_sample_seed_shifted_totals += short_seq.samples_counts
                short_seq.iso_5p = mod_thing.iso_5p
            else:
                short_seq.iso_5p = 0

            ## 3p Alt cut? / iso_3p?

            if 'is_3p_alt_cut' in disagreements:
                mod = 'is_3p_alt_cut'
                pmsg(disagreement_error_msg.format(mod, mod, main_seq,
                        short_seq, main_seq.genomic_locations,
                        short_seq.genomic_locations))
            elif mod_thing.is_3p_alt_cut:
                self.per_sample_3p_alt_cut_totals += short_seq.samples_counts
                short_seq.iso_3p = mod_thing.iso_3p
            else:
                short_seq.iso_3p = 0

            ## 3p Mismatch? / iso_add?

            if 'is_3p_mismatch' in disagreements:
                mod = 'is_3p_mismatch'
                pmsg(disagreement_error_msg.format(mod, mod, main_seq,
                        short_seq, main_seq.genomic_locations,
                        short_seq.genomic_locations))
            elif mod_thing.is_3p_mismatch:
                self.per_sample_3p_mismatch_totals += short_seq.samples_counts
                short_seq.iso_add = mod_thing.iso_add
            else:
                short_seq.iso_add = 0

            # Editing....
            # Note: Logically, since we only report iso_snp_* and is_*_edited for
            #       gen_loc_bins with a DESIGNATION_ONE binstarter, there should
            #       never be a disagreement in an edited field.  Nonetheless,
            #       for consistency and just in case of logic error, logging this.
            #       See counts listed in note above.

            ## Seed Edited?

            if 'is_seed_edited' in disagreements:
                mod = 'is_seed_edited'
                pmsg(disagreement_error_msg.format(mod, mod, main_seq,
                        short_seq, main_seq.genomic_locations,
                        short_seq.genomic_locations))
            elif mod_thing.is_seed_edited:
                self.per_sample_seed_edited_totals += short_seq.samples_counts

            ## iso_snp_seed?

            if 'iso_snp_seed' in disagreements:
                mod = 'iso_snp_seed'
                pmsg(disagreement_error_msg.format(mod, mod, main_seq,
                        short_seq, main_seq.genomic_locations,
                        short_seq.genomic_locations))
            elif mod_thing.iso_snp_seed:
                short_seq.iso_snp_seed = True
            else:
                short_seq.iso_snp_seed = False

            ## iso_snp_central_offset?

            if 'iso_snp_central_offset' in disagreements:
                mod = 'iso_snp_central_offset'
                pmsg(disagreement_error_msg.format(mod, mod, main_seq,
                        short_seq, main_seq.genomic_locations,
                        short_seq.genomic_locations))
            elif mod_thing.iso_snp_central_offset:
                short_seq.iso_snp_central_offset = True
            else:
                short_seq.iso_snp_central_offset = False

            ## Central Edited? / iso_snp_central?

            if 'is_central_edited' in disagreements:
                mod = 'is_central_edited'
                pmsg(disagreement_error_msg.format(mod, mod, main_seq,
                        short_seq, main_seq.genomic_locations,
                        short_seq.genomic_locations))
            elif mod_thing.is_central_edited:
                self.per_sample_central_edited_totals += short_seq.samples_counts
                short_seq.iso_snp_central = True
            else:
                short_seq.iso_snp_central = False

            ## 3p Supplementary Edited? / iso_central_supp?

            if 'is_3p_supplementary_edited' in disagreements:
                mod = 'is_3p_supplementary_edited'
                pmsg(disagreement_error_msg.format(mod, mod, main_seq,
                        short_seq, main_seq.genomic_locations,
                        short_seq.genomic_locations))
            elif mod_thing.is_3p_supplementary_edited:
                self.per_sample_3p_supplementary_edited_totals += short_seq.samples_counts
                short_seq.iso_central_supp = True
            else:
                short_seq.iso_central_supp = False

            ## Other Edited? / iso_snp?

            other_edited = False

            ## 3p first...

            if 'is_other_edited_3p' in disagreements:
                mod = 'is_other_edited_3p'
                pmsg(disagreement_error_msg.format(mod, mod, main_seq,
                        short_seq, main_seq.genomic_locations,
                        short_seq.genomic_locations))
            elif mod_thing.is_other_edited_3p:
                other_edited = True

            ## 5p now...

            if 'is_other_edited_5p' in disagreements:
                mod = 'is_other_edited_5p'
                pmsg(disagreement_error_msg.format(mod, mod, main_seq,
                        short_seq, main_seq.genomic_locations,
                        short_seq.genomic_locations))
            elif mod_thing.is_other_edited_5p:
                other_edited = True

            ## Finalize

            if other_edited:
                self.per_sample_other_edited_totals += short_seq.samples_counts
                short_seq.iso_snp = True
            else:
                short_seq.iso_snp = False


class ModificationThing(object):
    """What a terrible name.  Still, I kinda like it now.

    Feed two AlignmentExecutionHits to create a ModificationThing, one from a
    BinStarter, and one from a MemberSequence.  Then your newly created
    ModificationThing will tell you if the MemberSequence is:
    * seed shifted
    * mismatched on the 3p side (putatitve untemplated addition)
    * longer/shorter (put alt cut or post-cut trimmed)
    * seed edited
    * central edited
    * supplementary edited
    * "other" edited

    Update 20180321:
        Originally Prost had the 6 categories of post transcriptional modifications.
        Then we added the 8 miRTop categories, which mostly overlap with the
        original 6 Prost categories, but not quite.  Here's where we stand today:

            iso_5p  == is_seed_shifted
            iso_3p  == is_3p_alt_cut
            iso_add == is_3p_mismatch

            is_seed_edited =                editing in NTs 2-8
            is_central_edited =             editing in NTs 9-12
            is_3p_supplementary_edited =    editing in NTs 13-17
            is_other_edited_5p =            editing in NT 1
            is_other_edited_3p =            editing in NTs 18-lastMatchingNT

            iso_snp_seed =                  editing in NTs 2-7
            iso_snp_central_offset =        editing in NTs 8
            iso_snp_central =               editing in NTs 9-12
            iso_snp_supp =                  editing in NTs 13-17
            iso_snp =                       editing anywhere else
            iso_snp_5p =                    editing in NT 1
            iso_snp_3p =                    editing in NTs 18-lastMatchingNT

    Examples:

        Example1::

            1234567890123456789012345
            CCCCCTTTTTGGGGGAAAAACCCCC    main_hit : 25=
                      X    XX            mismatches wrt main_hit
         AAACCCCCTTTTTCGGGGTGAAACCCCCAA
         123456789012345678901234567890  member   : 1=1X10=1X4=2X10=
           X          X    XX            mismatches wrt genome

            mismatches_coords_ref_main   = (0,11,16,17)
            mismatches_coords_ref_member = (2,13,18,19)
            offset = -2

        Example2::

            1234567890123456789012345
            CCCCCTTTTTGGGGGAAAAACCCCC    main_hit : 25=
                      X    XX            mismatches wrt main_hit
              CCCTTTTTCGGGGTGAAACCC
              123456789012345678901      member   : 8=1X4=2X6=
                      X    XX            mismatches wrt genome

            mismatches_coords_ref_main   = (11,16,17)
            mismatches_coords_ref_member = ( 9,14,15)
            offset = 2

    """

    __slots__ = (   'is_seed_shifted',
                    'is_3p_alt_cut',
                    'is_3p_mismatch',
                    'is_seed_edited',
                    'is_central_edited',
                    'is_3p_supplementary_edited',
                    'is_other_edited_5p',
                    'is_other_edited_3p',
                    'iso_5p',
                    'iso_3p',
                    'iso_add',
                    'iso_snp_seed',
                    'iso_snp_central_offset',
                    'iso_snp_central',
                    'iso_snp_supp',
                    'iso_snp_5p',
                    'iso_snp_3p', )

    @classmethod
    def disagreements(self, modification_things):
        """Given a list of ModificationThings, determine if they disagree about
        any of their properties (and if they do, return which properties are in
        disagreement).

        Arguments:
            modification_things ([ModificationThing]): A list of
                ModificationThings that will be compared to one another to
                determine if there are any disagreements between their
                properties.

        Returns:
            (str,): A tuple of properties (as strings) on which
                modification_things disagree.

        """
        vals = {}
        disagreements_list = []
        for prop in self.__slots__:
            vals[prop] = set()
        for mod_thing in modification_things:
            for prop in self.__slots__:
                vals[prop].add(getattr(mod_thing, prop))
        for prop in self.__slots__:
            if len(vals[prop]) > 1:
                # Above is equivalent to:
                # if sorted(vals[prop]) == [False, True]:
                disagreements_list.append(prop)
        return tuple(sorted(disagreements_list))

        # O(n^2) version
        #disagreements_list = set()
        #for prop in self.__slots__:
        #    for i in xrange(len(modification_things)):
        #        for j in xrange(i+1, len(modification_things)):
        #            m1 = modification_things[i]
        #            m2 = modification_things[j]
        #            if getattr(m1, prop) != getattr(m2, prop):
        #                disagreements_list.add(prop)
        #return tuple(sorted(disagreements_list))

    def __init__(self, main_hit, member_hit, main_seq, member_seq):
        """
        Arguments:
            main_hit (SamExtendedCigarAlignmentExecutionHit): The main hit
                (i.e. a gen_loc from a BinStarter) which we are using as a
                proxy for the reference genome (since it is full_length 100%
                identity).
            member_hit (SamExtendedCigarAlignmentExecutionHit): The member hit
                (i.e. a gen_loc from a member of a GenLocBin) from which we are
                retrieving mismatches w/r/t the main_hit.
            main_seq (ShortSeq): The ShortSeq that owns the main_hit.
                Yes, a little crusty, but we're adding this in for miRTop at
                a later date.
            member_seq (ShortSeq): The ShortSeq that owns the member_hit.
                Yes, a little crusty, but we're adding this in for miRTop at
                a later date.
        """
        # Initialize:
        self.is_seed_shifted = False
        self.is_3p_alt_cut = False
        self.is_3p_mismatch = False
        self.is_seed_edited = False
        self.is_central_edited = False
        self.is_3p_supplementary_edited = False
        self.is_other_edited_5p = False
        self.is_other_edited_3p = False
        self.iso_5p = 0
        self.iso_3p = 0
        self.iso_add = 0
        self.iso_snp_seed = False
        self.iso_snp_central_offset = False
        self.iso_snp_central = False
        self.iso_snp_supp = False
        self.iso_snp_5p = False
        self.iso_snp_3p = False


        # Same linkage group, same strand
        assert(main_hit.lg == member_hit.lg)
        assert(main_hit.on_minus_strand == member_hit.on_minus_strand)
        # MainSeq is always a 1
        assert(main_hit.has_100_percent_core_identity)
        assert(main_hit.is_full_length)
        # No indels!
        assert(not main_hit.has_indelnts)
        assert(not member_hit.has_indelnts)

        mm_coords = self._mismatches_coords(main_hit, member_hit)

        # Now set it in stone and discard mm_coords:
        self._set_properties(mm_coords, main_hit, member_hit, main_seq, member_seq)

    def _set_properties(self, mm_coords, main_hit, member_hit, main_seq, member_seq):
        """Instead of creating methods for each property, we simply
        precalculate them here.

        This method calculates and caches the following properties:
            * is_seed_shifted - If a member hit shares a the same 5p start as
                the main hit, then it is *not* seed shifted.
            * is_3p_mismatch - If a member hit has mismatches on its 3prime
                end, then it has a 3p mismatch.  Yup.
            * is_3p_alt_cut - If the last matching 3p nucleotide in the member
                hit has the same position as the last matching 3p nucleotide in
                the main hit, then it is *not* 3p alternatively cut.
            * is_seed_edited - If there a missatch in the seed region
                it is seed edited.
            * is_central_edited - If there a missatch in the central region
                (NTs 9-12) it is central edited.
            * is_3p_supplementary_edited - If there a missatch in the 3 prime
                supplementary region (NTs 9-12) it is supplementary edited.
            * is_other_edited_5p - If there a missatch at NT 1, then it
                is other_edited_5p.
            * is_other_edited_3p - If there a missatch at NT 18 until the last
                matching NT, then it is other_edited_3p.
            * Plus all the iso_* fields

        Arguments:
            mm_coords (int,): A list of coordinates of mismatches found in the
                member_hit w/r/t the main_hit (i.e. the coordinates returned are
                indexed relative to the main_hit).
            main_hit (SamExtendedCigarAlignmentExecutionHit): The main hit
                (i.e. a gen_loc from a BinStarter) which we are using as a
                proxy for the reference genome (since it is full_length 100%
                identity).
            member_hit (SamExtendedCigarAlignmentExecutionHit): The member hit
                (i.e. a gen_loc from a member of a GenLocBin) from which we are
                retrieving mismatches w/r/t the main_hit.
            main_seq (ShortSeq): The ShortSeq that owns the main hit.
            member_seq (ShortSeq): The ShortSeq that owns the member hit.
            """

        ## Bail and log if the BinStarter aligns to a region with an 'N'
        if any(t[1] == 'M' for t in main_hit.cigar_tokens):
            raise ModificationThingEncounteredNAlignment

        ## Seed shifted?
        self.is_seed_shifted = (
                not (member_hit.reference_start_with_clips ==
                        main_hit.reference_start_with_clips))

        ## iso_5p ?
        if self.is_seed_shifted:
            self.iso_5p = (main_hit.reference_start_with_clips - \
                    member_hit.reference_start_with_clips)
            if member_hit.on_minus_strand:
                self.iso_5p = -self.iso_5p

        ## 3p mismatch?
        # TODO: consider adding 'S' and 'H' here as well?
        # Recall, cigar_tokens_5pto3p look like this: ((1, 'X'), (18, '='), (2, 'X'))
        if member_hit.cigar_tokens_5pto3p[-1][1] == 'X':
            self.is_3p_mismatch = True
            self.iso_add = member_hit.cigar_tokens_5pto3p[-1][0]
        elif member_hit.cigar_tokens_5pto3p[-1][1] == '=':
            pass
        elif member_hit.cigar_tokens_5pto3p[-1][1] == 'M':
            # Encountered an M in one of the alignments. Ignore and
            # log.
            raise ModificationThingEncounteredNAlignment
        else:
            raise MirModificationCalculationException, (
                "token_kind '{}' not allowed (E1)\n"
                "main_hit = {}\n"
                "member_hit = {}".format(
                    member_hit.cigar_tokens_5pto3p[-1][1], main_hit,
                    member_hit))

        ## These variables ('offset' and the 'member_last_matching_nt_3p_end')
        ## are needed for subsequent classifications (e.g. 3p alt cut and other_edited and so forth...)

        ## Get offset: How many nts on 5p side the member_hit is offset from main_hit.
        offset = (member_hit.reference_start_with_clips -
                main_hit.reference_start_with_clips)
        if main_hit.on_minus_strand:
            offset = -offset

        ## Last matching 3p nucleotide for main:
        ## (e.g. "22" for example)
        if main_hit.on_minus_strand:
            main_last_matching_nt_3p_end = (
                    main_hit.reference_start_with_clips
                    - main_hit.reference_end_with_clips
                    + 1)

        else:
            main_last_matching_nt_3p_end = (
                    main_hit.reference_end_with_clips
                    - main_hit.reference_start_with_clips
                    + 1)

        ## Last matching 3p nucleotide for member:
        ## (e.g. "22" for example)
        if member_hit.on_minus_strand:
            member_last_matching_nt_3p_end = (
                member_hit.reference_start_with_clips
                - member_hit.reference_end_with_clips
                + 1
                - member_hit.num_3p_mismatches
                + offset)
        else:
            member_last_matching_nt_3p_end = (
                member_hit.reference_end_with_clips
                - member_hit.reference_start_with_clips
                + 1
                - member_hit.num_3p_mismatches
                + offset)

        ## 3p alt cut?
        self.is_3p_alt_cut = not (
                main_last_matching_nt_3p_end == member_last_matching_nt_3p_end)

        ## iso_3p
        if self.is_3p_alt_cut:
            self.iso_3p = (member_last_matching_nt_3p_end -
                    main_last_matching_nt_3p_end)



        ## is_*_edited and iso_snp_* ?
        #
        # We note that since we only report is_*_edited iso_snp_* for sequences
        # that are a member of a gen_loc_bin with a DESIGNATION_ONE binstarter,
        # that there can never be a disagreement between mod_things w/r/t any
        # of the "edited" / "iso_snp_*" categories.
        #
        # iso_snp_seed           -> SNP at NTs 2-7
        # iso_snp_central_offset -> SNP at NT 8
        # iso_snp_central        -> SNP at NTs 9-12
        # iso_snp_supp           -> SNP at NTs 13-17
        # iso_snp                -> SNP anywhere else.
        # iso_snp_5p             -> SNP at NT 1
        # iso_snp_3p             -> SNP at NT 18 onwards to last matching NT
        #
        # Note also that we keep our definition of is_seed_edited to be NTs 2-8.

        for mm in mm_coords:
            if (mm >= 2 and mm <= 8):
                self.is_seed_edited = True
            if (mm >= 2 and mm <= 7):
                self.iso_snp_seed = True
            if (mm == 8):
                self.iso_snp_central_offset = True
            if (mm >= 9 and mm <= 12):
                self.is_central_edited = True
                self.iso_snp_central = True
            if (mm >= 13 and mm <= 17):
                self.is_3p_supplementary_edited = True
                self.iso_snp_supp = True
            if (mm == 1):
                self.is_other_edited_5p = True
                self.iso_snp_5p = True
            ## other edited on 3p end?
            # If the mismatch is located outside of the seed region (nts 2-8),
            # outside the central region (nts 9-12), or outside the
            # 3p-supplemenary-region (nts 13-17), then it is "other edited".
            # To distinguish is_other_edited from is_3p_mismatch, we force
            # is_other_edited to only consider nucleotides from nt17 until the
            # last matching nucleotide on the 3p end (basically, the last
            # matching nucleotide is the "border" between is_other_edited and
            # is_3p_mismatch).
            if (mm >= 18 and mm < member_last_matching_nt_3p_end):
                self.is_other_edited_3p = True
                self.iso_snp_3p = True

    @classmethod
    def _mismatches_coords(self, main_hit, member_hit):
        """Coordinates of mismatches w/r/t the main_hit.

        Arguments:
            main_hit (SamExtendedCigarAlignmentExecutionHit): The main hit
                (i.e. a gen_loc from a BinStarter) which we are using as a
                proxy for the reference genome (since it is full_length 100%
                identity).
            member_hit (SamExtendedCigarAlignmentExecutionHit): The member hit
                (i.e. a gen_loc from a member of a GenLocBin) from which we are
                retrieving mismatches w/r/t the main_hit.

        Returns:
            (int,): A list of coordinates of mismatches found in the member_hit
                w/r/t the main_hit (i.e. the coordinates returned are indexed
                relative to the main_hit).
        """

        curr_pos = 1
        mm_coords = []

        # How many nts on 5p side the member_hit is offset from main_hit.
        offset = (member_hit.reference_start_with_clips -
                main_hit.reference_start_with_clips)
        if main_hit.on_minus_strand:
            offset = -offset

        # Get cigar tokens (5p-to-3p) as mutable list of lists
        main_cigar_tokens = [list(ct) for ct in
                main_hit.cigar_tokens_5pto3p]
        member_cigar_tokens = [list(ct) for ct in
                member_hit.cigar_tokens_5pto3p]

        # Trim any 5p offset from appropriate hit's cigar_tokens.
        # (In Example1, member_hit leading 'AA' is trimmed, in Example2,
        # main_hit leading 'CC' is trimmed.)
        if offset > 0:
            # Trim main.
            tokens_trimmed = self._trim_left_tokens(main_cigar_tokens, offset)
            curr_pos = len(tokens_trimmed) + 1
        elif offset < 0:
            # Trim member.  Additionally, record any mismatches found in
            # member's 5p overhang.
            # NOTE:     - possibly allow 'S' or 'H' or 'P?? (no)' after
            #             thoughtful consideration.  For now, assert otherwise.
            tokens_trimmed = self._trim_left_tokens(
                    member_cigar_tokens, -offset)
            for i, token_kind in enumerate(tokens_trimmed):
                if token_kind == 'X':
                    mm_coords.append(i + offset + 1)
                elif token_kind == '=':
                    pass
                elif token_kind == 'M':
                    # Encountered an M in one of the alignments. Ignore and
                    # log.
                    raise ModificationThingEncounteredNAlignment
                else:
                    raise MirModificationCalculationException, (
                        "token_kind '{}' not allowed (E2)\n"
                        "main_hit = {}\n"
                        "member_hit = {}".format(
                            token_kind, main_hit, member_hit))

        # Walk down member_cigar_tokens, recording coord(s) of mismatch(es)
        for token in member_cigar_tokens:
            if token[1] == 'X':
                for i in xrange(token[0]):
                    mm_coords.append(curr_pos + i)
            elif token[1] == '=':
                pass
            elif token[1] == 'M':
                # Encountered an M in one of the alignments. Ignore and
                # log.
                raise ModificationThingEncounteredNAlignment
            else:
                raise MirModificationCalculationException, (
                    "token_kind '{}' not allowed (E3)\n"
                    "main_hit = {}\n"
                    "member_hit = {}".format(token[1], main_hit, member_hit))

            curr_pos += token[0]

        return tuple(mm_coords)

    @staticmethod
    def _trim_left_tokens(cigar_tokens, num_nts):
        """Trim num_nts tokens from the "left" side (equivalent to a perl array
        "shift").  If you pass in cigar_tokens that read 5p-to-3p, it will trim
        from the 5p side (and conversely will trim from 3p side for 3p-to-5p
        tokens).

        This method modifies the argument cigar_tokens.

        Examples:
            Example1:
                cigar_tokens = [[17, '='], [1, 'X'], [5, '=']]
                _trim_left_tokens(cigar_tokens, 2)
                # Now: cigar_tokens == [[15, '='], [1, 'X'], [5, '=']]
            Example2:
                cigar_tokens = [[1, '='], [1, 'X'], [17, '=']]
                _trim_left_tokens(cigar_tokens, 3)
                # Now cigar_tokens == [[16, '=']]

        Note:
            This would be inefficient if we were trimming large numbers of 5p
            nts (because it walks down the cigar tokens one nt at a time).
            However, this algorithm is probably more efficient for very small
            num_nts, which represents vast majority of cases.

        Arguments:
            cigar_tokens ([[int,String]]): The list version of cigar_tokens.  For
                example, [[17, '='], [1, 'X'], [5, '=']] represents the cigar
                string "17=1X5=".  See parse_cigar_tokens() for more info.
            num_nts (int): The number of nucleotides to trim off the
                cigar_tokens.
        Returns:
            String: A 'list' of tokens trimmed.  In Example1 above, it is "==",
                in Example2 is it "=X=".

        """
        assert(num_nts > 0)
        tokens_trimmed = ''
        for i in xrange(num_nts):
            if cigar_tokens[0][0] == 1:
                # First cigar_token is 1nt long, simply remove it
                token_trimmed = cigar_tokens.pop(0)
                tokens_trimmed += token_trimmed[1]
            else:
                # First cigar_token is > 1nt long, decrement:
                tokens_trimmed += cigar_tokens[0][1]
                cigar_tokens[0][0] -= 1
            assert(len(cigar_tokens) > 0)

        return tokens_trimmed


class Output(object):

    def _headers_annotations(self, conf, alignments,
            include_reverse_anno=False):
        """Create the headers for the annotations columns.

        Args:
            conf (Configuration): The Configuration singleton object.
            alignments (Alignments): The Alignments singleton.
            include_reverse_anno (bool): Whether or not to include the reverse
                annotations...

        Returns:
            [str]: A list of annotations headers.

        """
        header = []
        for annotation_alignment in alignments.annotation_alignments:
            annotation_cls = annotation_alignment.annotation_cls
            if annotation_cls == BiomartOtherRNAAnnotation:
                # Recall that BiomartOtherRNAAnnotations spans 2 cols.
                header.append('other_ncRNA')
                header.append('ncRNA_biotype')
            else:
                for kind in annotation_cls.Kind:
                    if kind.name == 'species':
                        if annotation_cls == MirbaseMirAnnotation:
                            header.append('{}_miRNA'.format(conf.general.species))
                        elif annotation_cls == MirbaseMirReverseAnnotation:
                            if include_reverse_anno:
                                header.append('{}_miRNA_rev'.format(conf.general.species))
                        elif annotation_cls == MirbaseHairpinAnnotation:
                            header.append('{}_hairpin'.format(conf.general.species))
                        else:
                            raise ControlFlowException, \
                                    "ERR415: Shouldn't be possible to reach here."
                    else:
                        if annotation_cls == MirbaseMirAnnotation:
                            header.append('other_species_miRNA')
                        elif annotation_cls == MirbaseMirReverseAnnotation:
                            if include_reverse_anno:
                                header.append('other_species_miRNA_rev')
                        elif annotation_cls == MirbaseHairpinAnnotation:
                            header.append('other_species_hairpin')
                        else:
                            raise ControlFlowException, \
                                    "ERR415: Shouldn't be possible to reach here."
        return header

    def tsv_isomirs(self, conf, samples, short_seqs, alignments, gen_loc_bins,
         seed_bins, annotation_bins):
        """Generates the uncompressed Prost output TSV file.

        Args:
            conf (Configuration): The Configuration singleton object.
            samples (Samples): The Samples singleton object.
            short_seqs (ShortSeqs): The ShortSeqs singleton object.
            alignments (Alignments): The Alignments singleton object.
            gen_loc_bins (GenLocBins): The GenLocBins singleton object.
            seed_bins (SeedBins): The SeedBins singleton object.
            annotation_bins (AnnotationBins): The AnnotationBins singleton object.

        """

        progress = Progress("Writing uncompressed isomiRs file", 100, len(short_seqs))
        with open(conf.general.output_file_isomirs, 'w') as f:

            ### Header ###

            header = ("Sequence Seed Loc_idx Seed_idx Anno_idx Locations "
                "CIGARs_5pto3p Designations BinStarter".split())

            # Extract indexes by which to sort.  We'll be sorting by binstarter,
            # designation, then sequence.
            idx_binstarter = header.index('BinStarter')
            idx_designation = header.index('Designations')
            idx_sequence = header.index('Sequence')

            # Samples
            for sample_name in samples.iterkeys():
                header.append(sample_name)
            for sample_name in samples.iterkeys():
                header.append("{}_norm".format(sample_name))

            # Annotations
            header += self._headers_annotations(conf, alignments, True)

            # miRTop columns
            header += "iso_5p iso_3p iso_add iso_snp_seed iso_snp_central_offset iso_snp_central iso_snp_supp iso_snp".split()

            # Write
            header = "{}\n".format("\t".join(header))
            f.write(header)

            ### Data ###

            rows = []

            # abbreviations for readability
            max_locs_to_report = conf.general.max_locations_to_report

            for short_seq in short_seqs.itervalues():
                progress.progress()

                if short_seq.is_no_hit:
                    continue

                # Get bins
                gen_loc_bin = gen_loc_bins.get_bin_by_seq(short_seq.seq_str)
                seed_bin = seed_bins.get_bin_by_seq(short_seq.seq_str)
                annotation_bin = annotation_bins.get_bin_by_seq(short_seq.seq_str)

                row = [short_seq.seq_str, short_seq.seed]

                # Bin indexes
                if gen_loc_bin is None:
                    assert (short_seq.ambiguous
                            or short_seq.is_no_hit
                            or short_seq.max_locations_allowed_exceeded)
                    row.append("")
                else:
                    row.append(gen_loc_bin.idx)
                if seed_bin is None:
                    assert(short_seq.is_no_hit or annotation_bin is None)
                    row.append("")
                else:
                    row.append(seed_bin.idx)
                if annotation_bin is None:
                    row.append("")
                else:
                    row.append(annotation_bin.idx)

                # Locations
                gen_locs = short_seq.genomic_locations
                count_gen_locs = short_seq.count_genomic_locations
                if short_seq.max_locations_allowed_exceeded:
                    row.append("MLAE:{}".format(count_gen_locs))
                elif count_gen_locs > max_locs_to_report:
                    row.append("TML:{}".format(count_gen_locs))
                else:
                    locs = []
                    for loc in gen_locs:
                        locs.append("{}:{}-{}".format(loc.lg, loc.start, loc.end))
                    row.append(";".join(locs))

                # Cigar strings
                if short_seq.max_locations_allowed_exceeded:
                    row.append("MLAE:{}".format(count_gen_locs))
                elif count_gen_locs > max_locs_to_report:
                    row.append("TML:{}".format(count_gen_locs))
                else:
                    locs = []
                    for loc in gen_locs:
                        locs.append("{}".format(loc.cigar5pto3p))
                    row.append(";".join(locs))

                # Designations
                if short_seq.max_locations_allowed_exceeded:
                    # report 'MLAE'
                    desig = "MLAE"
                elif ((count_gen_locs > max_locs_to_report)
                        or short_seq.designation_integer < 3
                        or short_seq.is_no_hit):
                    # report just the designation int for 1's, 2's, and TMLs
                    desig = str(short_seq.designation_integer)
                else:
                    # report all the designations for >= 3's (if not too many)
                    desig = ";".join((l.designation for l in gen_locs))
                row.append(desig)

                # Bin starter
                if short_seq.ambiguous:
                    row.append('ambiguous: {}'.format(
                        ','.join(short_seq.ambiguous)))
                elif (short_seq.is_no_hit
                        or short_seq.max_locations_allowed_exceeded):
                    row.append("")
                else:
                    row.append(gen_loc_bin.main_short_seq_str)

                # Samples
                for c in short_seq.samples_counts.counts():
                    row.append(str(c))
                for n in short_seq.samples_counts.norms():
                    row.append(str(n))

                # annotations
                # We walk through the sequence's annotations list for every
                # (annotation_class, kind) pair and collect them.  The
                # exception is BiomartOtherRNAAnnotation, for which we create
                # two cells (across two columns) instead of one.
                for annotation_alignment in alignments.annotation_alignments:
                    annotation_cls = annotation_alignment.annotation_cls
                    if annotation_cls == BiomartOtherRNAAnnotation:
                        annos = annotations_of(BiomartOtherRNAAnnotation,
                            short_seqs, None, short_seq.seq_str)
                        names = []
                        biotypes = []
                        for anno in sorted(annos, key=operator.attrgetter('name')):
                            if anno.name not in names:
                                names.append(anno.name)
                                biotypes.append(anno.biotype)
                        row.append(",".join(names))
                        row.append(",".join(biotypes))
                    else:
                        for kind in annotation_cls.Kind:
                            annos = []
                            for anno in short_seq.annotations:
                                if (anno.__class__ == annotation_cls
                                        and anno.kind == kind):
                                    annos.append(anno.name)
                            row.append(",".join(sorted(annos)))

                # miRTop columns
                mirtops = [short_seq.iso_5p, short_seq.iso_3p,
                        short_seq.iso_add, short_seq.iso_snp_seed,
                        short_seq.iso_snp_central_offset,
                        short_seq.iso_snp_central, short_seq.iso_snp_supp,
                        short_seq.iso_snp]

                mirtops = [str(m) for m in mirtops]
                row += mirtops

                # Save
                rows.append(row)

            # Sort: by binstarter, designation, then sequence
            rows = sorted(rows, key=operator.itemgetter(
                    idx_binstarter, idx_designation, idx_sequence))

            # Write
            for row in rows:
                row = "{}\n".format("\t".join(row))
                f.write(row)

        progress.done()

    def tsv_no_hits(self, conf, samples, short_seqs, alignments):
        """Generates the uncompressed Prost output TSV file of sequences that do
        not hit to the genome.

        Args:
            conf (Configuration): The Configuration singleton object.
            samples (Samples): The Samples singleton object.
            short_seqs (ShortSeqs): The ShortSeqs singleton object.
            alignments (Alignments): The Alignments singleton object.

        """

        progress = Progress("Writing no_hits output file", 100, len(short_seqs))
        with open(conf.general.output_file_no_hits, 'w') as f:

            ### Header ###

            header = "Sequence Seed".split()

            # Samples
            sample_names = []
            for sample_name in samples.iterkeys():
                sample_names.append(sample_name)
            header += sample_names

            # Annotations
            header += self._headers_annotations(conf, alignments, True)

            # Include post-filtered locations/cigars
            header += "Locations CIGARs_5pto3p soft_clipped?".split()

            # Extract indexes by which to sort.  We'll be sorting first by
            # annotations, then by the sum of unnormalized counts (higest
            # expression goes first).
            anno_idxs = [i for i, h in enumerate(header) if
                    re.search('(_miRNA|_hairpin|_ncRNA)$', h)]
            idx_anno_first, idx_anno_last = anno_idxs[0], anno_idxs[-1]
            expr_idxs = [i for i, h in enumerate(header) if h in sample_names]
            idx_expr_first, idx_expr_last = expr_idxs[0], expr_idxs[-1]

            # Write
            header = "{}\n".format("\t".join(header))
            f.write(header)

            ### Data ###

            rows = []

            # abbreviations for readability
            max_locs_allowed = conf.general.max_locations_allowed
            max_locs_to_report = conf.general.max_locations_to_report

            for short_seq in short_seqs.itervalues():
                progress.progress()

                if not short_seq.is_no_hit:
                    continue

                row = [short_seq.seq_str, short_seq.seed]

                # Designations
                desig = str(short_seq.designation_integer)

                # Samples
                row += short_seq.samples_counts.counts()

                # annotations
                # We walk through the sequence's annotations list for every
                # (annotation_class, kind) pair and collect them.  The
                # exception is BiomartOtherRNAAnnotation, for which we create
                # two cells (across two columns) instead of one.
                for annotation_alignment in alignments.annotation_alignments:
                    annotation_cls = annotation_alignment.annotation_cls
                    if annotation_cls == BiomartOtherRNAAnnotation:
                        annos = []
                        biotypes = []
                        for anno in short_seq.annotations:
                            if (anno.__class__ == annotation_cls):
                                annos.append(anno.name)
                                biotypes.append(anno.biotype)
                        row.append(",".join(annos))
                        row.append(",".join(biotypes))
                    else:
                        for kind in annotation_cls.Kind:
                            annos = []
                            for anno in short_seq.annotations:
                                if (anno.__class__ == annotation_cls
                                        and anno.kind == kind):
                                    annos.append(anno.name)
                            row.append(",".join(annos))


                # Prost-post-filtered-out Locations
                no_hit_gen_locs = short_seq.no_hit_genomic_locations
                count_no_hit_gen_locs = short_seq.count_no_hit_genomic_locations

                if count_no_hit_gen_locs > max_locs_allowed:
                    row.append("MLAE:{}".format(count_no_hit_gen_locs))
                elif count_no_hit_gen_locs > max_locs_to_report:
                    row.append("TML:{}".format(count_no_hit_gen_locs))
                elif count_no_hit_gen_locs == 0:
                    # no "no_hit" locations...
                    row.append("")
                else:
                    locs = []
                    for loc in no_hit_gen_locs:
                        locs.append("{}:{}-{}".format(loc.lg, loc.start, loc.end))
                    row.append(";".join(locs))

                # Cigar strings for those locations
                if count_no_hit_gen_locs > max_locs_allowed:
                    row.append("MLAE:{}".format(count_no_hit_gen_locs))
                elif count_no_hit_gen_locs > max_locs_to_report:
                    row.append("TML:{}".format(count_no_hit_gen_locs))
                elif count_no_hit_gen_locs == 0:
                    # no "no_hit" locations...
                    row.append("")
                else:
                    locs = []
                    for loc in no_hit_gen_locs:
                        locs.append("{}".format(loc.cigar5pto3p))
                    row.append(";".join(locs))

                # Soft Clipped?
                row.append(loc.is_soft_clipped)

                # Save
                rows.append(row)

            # Sort: first by annotations, then by sum of unnormalized counts
            rows = sorted(rows, key=lambda row: (
                    sum(row[idx_expr_first:idx_expr_last+1])), reverse=True)
            rows = sorted(rows, key=lambda row: (
                    row[idx_anno_first:idx_anno_last+1]), reverse=True)

            # Write
            for row in rows:
                row = "{}\n".format("\t".join(str(e) for e in row))
                f.write(row)

        progress.done()

    def tsv_comp_by_gen_loc(self, conf, samples, short_seqs, bins, alignments):
        """Generates the compressed by GenomicLocation Prost output TSV file.

        Args:
            conf (Configuration): The Configuration singleton object.
            samples (Samples): The Samples singleton object.
            short_seqs (ShortSeqs): The ShortSeqs singleton object.
            bins (GenLocBins): The GenLocBins singleton object.
            alignments (Alignments): The Alignments singleton object.
        """

        progress = Progress("Writing compressed output file", 100, len(bins))
        with open(conf.general.output_file_comp_by_gen_loc, 'w') as f:

            ### Header ###

            header = "Loc_idx BinStarter Locations CIGARs_5pto3p Designations".split()

            # Annotations
            header += self._headers_annotations(conf, alignments, True)

            # "Extract" indexes by which to sort.  We'll be sorting first by the
            # the designation (we want "1s" on top), then in_species_hairpin,
            # followed by the reverse of the in_species_miRNA, then
            # other_species_hairpin, then other_species_miRNA (reverse). The
            # reason for this sorting is so that roughly the 5p will come before
            # the 3p (it's not perfect, but it's good enough).
            #
            # Note that this is HARDCODED because it's easier for this
            # particular case (the header names are not set in stone for now...)
            idx_designation = header.index('Designations')
            idx_anno_in_species_hairpin = 9
            idx_anno_in_species_miR = 5
            idx_anno_other_species_hairpin = 10
            idx_anno_other_species_miR = 6
            idx_anno_other_species_ncRNA = 11

            # Ambiguous this_species mir and hairpin annotations
            header.append('{}_miRNA_ambiguous?'.format(conf.general.species))
            header.append('{}_hairpin_ambiguous?'.format(conf.general.species))

            # Samples
            # i.e. sample, _norm, seed shifted, seed/3'suppl/other edited,
            # 3' alt_cut/mismatch
            for sample_name in samples.iterkeys():
                header.append(sample_name)
            for sample_name in samples.iterkeys():
                header.append("{}_norm".format(sample_name))
            for sample_name in samples.iterkeys():
                header.append("{}_%seed_shifted".format(sample_name))
            for sample_name in samples.iterkeys():
                header.append("{}_%seed_edited".format(sample_name))
            for sample_name in samples.iterkeys():
                header.append("{}_%central_edited".format(sample_name))
            for sample_name in samples.iterkeys():
                header.append("{}_%3'-supplementary_edited".format(sample_name))
            for sample_name in samples.iterkeys():
                header.append("{}_%3'-alternatively_cut".format(sample_name))
            for sample_name in samples.iterkeys():
                header.append("{}_%3'-mismatch".format(sample_name))
            for sample_name in samples.iterkeys():
                header.append("{}_%other_edited".format(sample_name))

            # Gapped
            header.append('Gapped?')

            # Write header
            header = "{}\n".format("\t".join(header))
            f.write(header)

            ### Data ###

            rows = []
            for bn in bins.all_bins_sorted_by_starting_seq():
                progress.progress()

                # Bin starting sequence
                starter_seq = short_seqs[bn.main_short_seq_str]

                # Add the bin starting sequence and bin index
                row = [bn.idx, starter_seq.seq_str]

                # Locations
                gen_locs = starter_seq.genomic_locations
                if len(gen_locs) > conf.general.max_locations_to_report:
                    row.append("TML:{}".format(len(gen_locs)))
                else:
                    locs = []
                    for loc in gen_locs:
                        locs.append("{}:{}-{}".format(loc.lg, loc.start, loc.end))
                    row.append(";".join(locs))

                # CIGAR strings
                if len(gen_locs) > conf.general.max_locations_to_report:
                    row.append("TML:{}".format(len(gen_locs)))
                else:
                    locs = []
                    for loc in gen_locs:
                        locs.append("{}".format(loc.cigar5pto3p))
                    row.append(";".join(locs))

                # Designations
                if (len(gen_locs) >
                        conf.general.max_locations_to_report
                        or starter_seq.designation_integer < 3):
                    # report just the designation int for 1's, 2's, and TMLs
                    desig = str(starter_seq.designation_integer)
                else:
                    # report all the designations for >= 3's (if not too many)
                    desig = ";".join((l.designation for l in gen_locs))
                row.append(desig)

                # Annotations
                #
                # We walk through every (annotation_class, kind) pair to create
                # the bin's annotations. (The exception is
                # BiomartOtherRNAAnnotation, for which we create two cells
                # instead of one (that is, BiomartOtherRNAAnnotation spans two
                # columns).
                #
                # See GenLocBin._mirbase_annotations for the specific algorithms
                # used to generate these lists of mirbase annotations.

                for annotation_alignment in alignments.annotation_alignments:
                    annotation_cls = annotation_alignment.annotation_cls
                    if annotation_cls == BiomartOtherRNAAnnotation:
                        names, biotypes = bn.biomart_other_rna_annotations(short_seqs)
                        row.append(",".join(names))
                        row.append(",".join(biotypes))
                    elif annotation_cls == MirbaseMirAnnotation:
                        for annos in bn.mirbase_mir_annotations(short_seqs):
                            row.append(",".join(a.name for a in annos))
                    elif annotation_cls == MirbaseMirReverseAnnotation:
                        for annos in bn.mirbase_mir_reverse_annotations(short_seqs):
                            row.append(",".join(a.name for a in annos))
                    elif annotation_cls == MirbaseHairpinAnnotation:
                        for annos in bn.mirbase_hairpin_annotations(short_seqs):
                            row.append(",".join(a.name for a in annos))

                # Ambiguous this_species mir and hairpin annotations?
                row.append(str(bn.mir_ambig(short_seqs)))
                row.append(str(bn.hairpin_ambig(short_seqs)))

                # Samples
                # i.e. sample, _norm, seed shifted, seed/3'suppl/other edited,
                # 3' alt_cut/mismatch
                row += [str(i) for i in bn.per_sample_read_totals.counts()]
                row += [str(i) for i in bn.per_sample_read_totals.norms()]
                na = (starter_seq.designation_integer >= DESIGNATION_THREE)
                row += [str(i) for i in bn.per_sample_perc_seed_shifted(na)]
                row += [str(i) for i in bn.per_sample_perc_seed_edited(na)]
                row += [str(i) for i in bn.per_sample_perc_central_edited(na)]
                row += [str(i) for i in bn.per_sample_perc_3p_supplementary_edited(na)]
                row += [str(i) for i in bn.per_sample_perc_3p_alt_cut(na)]
                row += [str(i) for i in bn.per_sample_perc_3p_mismatch(na)]
                row += [str(i) for i in bn.per_sample_perc_other_edited(na)]

                # Gapped
                row.append(str(bn.has_indelnts(short_seqs)))

                # Save
                rows.append(row)

            # Sort: In the following order:
            #   * designation (first character only)
            #   * in_species_hairpin
            #   * in_species_miRNA (reverse)
            #   * other_species_hairpin
            #   * other_species_miRNA (reverse)
            #   * other_species_ncRNA
            #   * designation (the entire string)
            # Use stability of sort (i.e. sort several times).

            # An initial sort by designation to get them in order.
            rows = sorted(rows, key=lambda row:(row[idx_designation]))
            idx = idx_anno_other_species_ncRNA
            rows = sorted(rows, key=lambda row: (
                    (row[idx] == "", row[idx]) ))
            idx = idx_anno_other_species_miR
            rows = sorted(rows, key=lambda row: (
                    (row[idx] == "", row[idx]) ), reverse=True)
            idx = idx_anno_other_species_hairpin
            rows = sorted(rows, key=lambda row: (
                    (row[idx] == "", row[idx]) ))
            idx = idx_anno_in_species_miR
            rows = sorted(rows, key=lambda row: (
                    (row[idx] == "", row[idx]) ), reverse=True)
            idx = idx_anno_in_species_hairpin
            rows = sorted(rows, key=lambda row: (
                    (row[idx] == "", row[idx]) ))
            rows = sorted(rows, key=lambda row:(row[idx_designation][0]))

            # Write
            for row in rows:
                row = "{}\n".format("\t".join(row))
                f.write(row)

        progress.done()

    def tsv_comp_by_seed(self, conf, samples, short_seqs, seed_bins,
                      gen_loc_bins, alignments):
        """Generates the seed-compressed Prost output TSV file.

        Args:
            conf (Configuration): The Configuration singleton object.
            samples (Samples): The Samples singleton object.
            short_seqs (ShortSeqs): The ShortSeqs singleton object.
            seed_bins (SeedBins): The SeedBins singleton object.
            gen_loc_bins (GenLocBins): The GenLocBins singleton object.
            alignments (Alignments): The Alignments singleton object.
        """

        progress = Progress("Writing seed-compressed output file", 100, len(seed_bins))
        with open(conf.general.output_file_comp_by_seed, 'w') as f:

            ### Header ###

            header = "Seed_idx Seed MainSequence".split()

            # Extract indexes by which to sort.  We'll be sorting by seed.
            idx_seed = header.index('Seed')

            # Annotations
            header += self._headers_annotations(conf, alignments)

            # Samples
            # i.e. sample, _norm
            for sample_name in samples.iterkeys():
                header.append(sample_name)
            for sample_name in samples.iterkeys():
                header.append("{}_norm".format(sample_name))

            # Gapped
            header.append('Gapped?')

            # Write header
            header = "{}\n".format("\t".join(header))
            f.write(header)

            ### Data ###

            rows = []
            for bn in seed_bins.all_bins_sorted_by_starting_seq():
                progress.progress()

                main_short_seq = short_seqs[bn.main_short_seq_str]

                # Bin seed
                seed = main_short_seq.seed

                # Add the seed and seed index
                row = [bn.idx, seed, bn.main_short_seq_str]

                # Annotations
                for annotation_alignment in alignments.annotation_alignments:
                    annotation_cls = annotation_alignment.annotation_cls
                    if annotation_cls == BiomartOtherRNAAnnotation:
                        names, biotypes = bn.biomart_other_rna_annotations(short_seqs)
                        row.append(",".join(names))
                        row.append(",".join(biotypes))
                    elif annotation_cls == MirbaseMirAnnotation:
                        for annos in bn.mirbase_mir_annotations(short_seqs,
                                                gen_loc_bins):
                            row.append(",".join(sorted(a.name for a in annos)))
                    elif annotation_cls == MirbaseHairpinAnnotation:
                        for annos in bn.mirbase_hairpin_annotations(short_seqs,
                                                gen_loc_bins):
                            row.append(",".join(sorted(a.name for a in annos)))

                # Samples
                # i.e. sample, _norm
                row += [str(i) for i in bn.per_sample_read_totals.counts()]
                row += [str(i) for i in bn.per_sample_read_totals.norms()]

                # Gapped
                row.append(str(bn.has_indelnts(short_seqs)))

                # Save
                rows.append(row)

            # Sort: by seed
            rows = sorted(rows, key=operator.itemgetter(idx_seed))

            # Write
            for row in rows:
                row = "{}\n".format("\t".join(row))
                f.write(row)

        progress.done()

    def tsv_comp_by_annotation(self, conf, samples, short_seqs, annotation_bins,
                            alignments):
        """Generates the annotation-compressed Prost output TSV file.

        Args:
            conf (Configuration): The Configuration singleton object.
            samples (Samples): The Samples singleton object.
            short_seqs (ShortSeqs): The ShortSeqs singleton object.
            annotation_bins (AnnotationBins): The AnnotationBins singleton object.
            alignments (Alignments): The Alignments singleton object.
        """

        progress = Progress("Writing annotation-compressed output file", 100,
                            len(annotation_bins))

        # Optional: Read in mature miR anno file, make a dict
        if conf.general.mature_mir_annotation_fasta:
            mature_mir_anno_dict = {}
            with open(conf.general.mature_mir_annotation_fasta, 'r') as f:
                while(True):
                    desc = f.readline().rstrip()
                    seq = f.readline().rstrip()
                    if not desc:
                        break
                    desc = desc.split()[0][1:]
                    mature_mir_anno_dict[desc] = seq

        with open(conf.general.output_file_comp_by_annotation, 'w') as f:

            ### Header ###
            header = "Anno_idx MainSequence".split()

            # MainSeqMatchesAnnotationFile
            if conf.general.mature_mir_annotation_fasta:
                header.append('MainSeqMatchesAnnotationFile')

            # Annotations
            header += self._headers_annotations(conf, alignments, True)

            # "Extract" indexes by which to sort.  We'll be sorting first by the
            # in_species_hairpin, followed by the reverse of the in_species_miRNA.
            # The reason for this sorting is so that roughly the 5p will come
            # before the 3p (it's not perfect, but it's good enough).
            #
            # Note that this is HARDCODED because it's easier for this
            # particular case (the header names are not set in stone for now...)
            idx_anno_in_species_hairpin = 5
            idx_anno_in_species_miR = 3

            # Samples
            # i.e. sample, _norm
            for sample_name in samples.iterkeys():
                header.append(sample_name)
            for sample_name in samples.iterkeys():
                header.append("{}_norm".format(sample_name))

            # Gapped
            header.append('Gapped?')

            # Write header
            header = "{}\n".format("\t".join(header))
            f.write(header)

            ### Data ###
            rows = []

            for bn in annotation_bins.all_bins_sorted_by_starting_seq():
                progress.progress()

                main_short_seq = short_seqs[bn.main_short_seq_str]

                # Add the annotations and bin index
                row = [bn.idx, bn.main_short_seq_str]

                # MainSeqMatchesAnnotationFile
                mir_species_annos = bn.mirbase_mir_annotations(short_seqs)[SPECIES_IDX]
                mir_species_rev_annos = bn.mirbase_mir_reverse_annotations(short_seqs)[SPECIES_IDX]
                if conf.general.mature_mir_annotation_fasta:
                    if len(mir_species_annos) > 1:
                        # multiple annotations
                        row.append('Multiple')
                    elif len(mir_species_annos) == 1:
                        anno = mir_species_annos[0].name
                        if mature_mir_anno_dict[anno] == bn.main_short_seq_str:
                            row.append('Yes')
                        else:
                            row.append('No')
                    elif len(mir_species_annos) == 0:
                        if len(mir_species_rev_annos) > 0:
                            row.append('Reverse')
                        else:
                            # none - should be impossible?
                            raise ControlFlowException, \
                                "ERR408: Shouldn't be possible to reach here."
                    else:
                        raise ControlFlowException, \
                            "ERR407: Shouldn't be possible to reach here."

                # Annotations
                for annotation_alignment in alignments.annotation_alignments:
                    annotation_cls = annotation_alignment.annotation_cls
                    if annotation_cls == BiomartOtherRNAAnnotation:
                        names, biotypes = bn.biomart_other_rna_annotations(short_seqs)
                        row.append(",".join(names))
                        row.append(",".join(biotypes))
                    elif annotation_cls == MirbaseMirAnnotation:
                        for annos in bn.mirbase_mir_annotations(short_seqs):
                            row.append(",".join(a.name for a in annos))
                    elif annotation_cls == MirbaseMirReverseAnnotation:
                        for annos in bn.mirbase_mir_reverse_annotations(short_seqs):
                            row.append(",".join(a.name for a in annos))
                    elif annotation_cls == MirbaseHairpinAnnotation:
                        for annos in bn.mirbase_hairpin_annotations(short_seqs):
                            row.append(",".join(a.name for a in annos))

                # Samples
                # i.e. sample, _norm
                row += [str(i) for i in bn.per_sample_read_totals.counts()]
                row += [str(i) for i in bn.per_sample_read_totals.norms()]

                # Gapped
                row.append(str(bn.has_indelnts(short_seqs)))

                # Save
                rows.append(row)

            # Sort: first by in_species_hairpin, then by in_species_miR
            # (reverse).  Use stability of sort (i.e. sort several times).
            rows = sorted(rows, key=lambda row: (
                    row[idx_anno_in_species_miR]), reverse=True)
            rows = sorted(rows, key=lambda row: (
                    row[idx_anno_in_species_hairpin]))

            # Write
            for row in rows:
                row = "{}\n".format("\t".join(row))
                f.write(row)

        progress.done()

    def tsv_mirror_mirs(self, conf, short_seqs, bins, alignments):
        """Generates the candidate mirror miRs Prost output TSV file.

        Args:
            conf (Configuration): The Configuration singleton object.
            samples (Samples): The Samples singleton object.
            short_seqs (ShortSeqs): The ShortSeqs singleton object.
            bins (GenLocBins): The GenLocBins singleton object.
            alignments (Alignments): The Alignments singleton object.
        """

        progress = Progress("Writing mirror miRNAs output file")
        with open(conf.general.output_file_mirror_mirs, 'w') as f:

            ### Header ###

            header =    "Loc_idx1 BinStarter1 Location1 " \
                        "Loc_idx2 BinStarter2 Location2".split()

            # Annotations
            anno_headers = self._headers_annotations(conf, alignments)
            anno_headers1 = []
            anno_headers2 = []
            for ah in anno_headers:
                anno_headers1.append("{}1".format(ah))
            for ah in anno_headers:
                anno_headers2.append("{}2".format(ah))
            header += anno_headers1
            header += anno_headers2

            # Write header
            header = "{}\n".format("\t".join(header))
            f.write(header)

            ### Data ###

            for (bn1, gen_loc1, bn2, gen_loc2) in bins.mirror_pairs:
                row = []

                # Bin starting sequences
                starter_seq1 = short_seqs[bn1.main_short_seq_str]
                starter_seq2 = short_seqs[bn2.main_short_seq_str]

                # Locations
                loc1 = "{}:{}-{}".format(
                                    gen_loc1.lg, gen_loc1.start, gen_loc1.end)
                loc2 = "{}:{}-{}".format(
                                    gen_loc2.lg, gen_loc2.start, gen_loc2.end)

                # Add the bin starting sequence, bin index, and potential
                # mirror_mir location.
                row += [bn1.idx, starter_seq1.seq_str, loc1]
                row += [bn2.idx, starter_seq2.seq_str, loc2]

                # Annotations
                #
                # We walk through every (annotation_class, kind) pair to create
                # the bin's annotations. (The exception is
                # BiomartOtherRNAAnnotation, for which we create two cells
                # instead of one (that is, BiomartOtherRNAAnnotation spans two
                # columns).
                #
                # See GenLocBin._mirbase_annotations for the specific algorithms
                # used to generate these lists of mirbase annotations.

                for bn in (bn1, bn2):
                    for annotation_alignment in alignments.annotation_alignments:
                        annotation_cls = annotation_alignment.annotation_cls
                        if annotation_cls == BiomartOtherRNAAnnotation:
                            names, biotypes = bn.biomart_other_rna_annotations(short_seqs)
                            row.append(",".join(names))
                            row.append(",".join(biotypes))
                        elif annotation_cls == MirbaseMirAnnotation:
                            for annos in bn.mirbase_mir_annotations(short_seqs):
                                row.append(",".join(a.name for a in annos))
                        elif annotation_cls == MirbaseHairpinAnnotation:
                            for annos in bn.mirbase_hairpin_annotations(short_seqs):
                                row.append(",".join(a.name for a in annos))

                # write
                row = "{}\n".format("\t".join(row))
                f.write(row)

        progress.done()

    def tsv_arm_switch(self, conf, samples, short_seqs, bins, alignments):
        """Generates the candidate arm switching Prost output TSV file.

        Arguments:
            conf (Configuration): The Configuration singleton object.
            samples (Samples): The Samples singleton object.
            short_seqs (ShortSeqs): The ShortSeqs singleton object.
            bins (AnnotationBins): The AnnotationBins singleton object.
            alignments (Alignments): The Alignments singleton object.

        """
        progress = Progress("Writing arm switch output file")
        with open(conf.general.output_file_arm_switch, 'w') as f:

            ### Header ###

            # Annotations mashup, prep
            # Basically, species_hairpin will be extracted and used untouched,
            # species_miR is redundant in this context and hence is deleted
            # altogether, while the remainder (i.e. other_(hairpin|miR),
            # other_ncRNA, and ncRNA_biotype) will have "_5p" and "_3p" appended
            # to each.
            anno_headers = self._headers_annotations(conf, alignments)
            anno_header_species_hairpin = anno_headers.pop(2)
            anno_header_species_miR = anno_headers.pop(0)
            anno_headers_5p = []
            anno_headers_3p = []
            for ah in anno_headers:
                anno_headers_5p.append("{}_5p".format(ah))
            for ah in anno_headers:
                anno_headers_3p.append("{}_3p".format(ah))

            # Initial Header
            header = """Anno_idx_5p Anno_idx_3p {} MainSequence_5p
                MainSequence_3p""".format(anno_header_species_hairpin).split()

            # Log fold changes
            # i.e. $(sample_name)_log_fold_change_expr_5p_3p
            for sample_name in samples.iterkeys():
                header.append("{}_log_fold_change_expr_5p_3p".format(sample_name))

            # Normalized Counts
            for sample_name in samples.iterkeys():
                header.append("{}_norm_5p".format(sample_name))
            for sample_name in samples.iterkeys():
                header.append("{}_norm_3p".format(sample_name))

            # Annotations
            header += anno_headers_5p
            header += anno_headers_3p

            # Extract indexes by which to sort.  We'll be sorting by the sum of
            # normalized expression (highest overall expression goes first):
            norm_idxs = \
                [i for i, h in enumerate(header) if re.search('_norm_(5|3)p$', h)]
            idx_norm_first, idx_norm_last = norm_idxs[0], norm_idxs[-1]

            # Write header
            header = "{}\n".format("\t".join(header))
            f.write(header)

            ### Data ###

            rows = []
            for (species_hairpin_anno_name, bn_5p, bn_3p) in bins.arm_5p_3p_pairs:
                row = []

                # Bin starting sequences
                starter_seq_5p = short_seqs[bn_5p.main_short_seq_str]
                starter_seq_3p = short_seqs[bn_3p.main_short_seq_str]

                # Add the annotation indexes, the species_hairpin name, and the
                # starting sequences for each bin
                row += [bn_5p.idx, bn_3p.idx, species_hairpin_anno_name,
                    starter_seq_5p.seq_str, starter_seq_3p.seq_str ]

                # Log Fold Changes (5p/3p)
                log_fold_changes = log_fold_changes_expr_5p_3p(
                        bn_5p.per_sample_read_totals,
                        bn_3p.per_sample_read_totals
                )
                row += log_fold_changes

                # Normalized Counts
                row += bn_5p.per_sample_read_totals.norms()
                row += bn_3p.per_sample_read_totals.norms()

                # Annotations
                #
                # We walk through every (annotation_class, kind) pair to create
                # the bin's annotations. (The exception is
                # BiomartOtherRNAAnnotation, for which we create two cells
                # instead of one (that is, BiomartOtherRNAAnnotation spans two
                # columns).
                #
                # See GenLocBin._mirbase_annotations for the specific algorithms
                # used to generate these lists of mirbase annotations.
                #
                # The only difference here is that we're ignoring in_species
                # mature and hairpin annotations.

                for bn in (bn_5p, bn_3p):
                    for annotation_alignment in alignments.annotation_alignments:
                        annotation_cls = annotation_alignment.annotation_cls
                        if annotation_cls == BiomartOtherRNAAnnotation:
                            names, biotypes = bn.biomart_other_rna_annotations(short_seqs)
                            row.append(",".join(names))
                            row.append(",".join(biotypes))
                        elif annotation_cls == MirbaseMirAnnotation:
                            # Skip species annos, only do other annos
                            o_idx = annotation_cls.Kind.other_species.value
                            annos = bn.mirbase_mir_annotations(short_seqs)[o_idx]
                            row.append(",".join(a.name for a in annos))
                        elif annotation_cls == MirbaseHairpinAnnotation:
                            # Skip species annos, only do other annos
                            o_idx = annotation_cls.Kind.other_species.value
                            annos = bn.mirbase_hairpin_annotations(short_seqs)[o_idx]
                            row.append(",".join(a.name for a in annos))

                # Save
                rows.append(row)

            # Sort: by sum of normalized counts
            rows = sorted(rows, key=lambda row:
                    (sum(row[idx_norm_first:idx_norm_last+1])), reverse=True)

            # Write
            for row in rows:
                row = "{}\n".format("\t".join(str(e) for e in row))
                f.write(row)

        progress.done()


#######################
### Debug functions ###
#######################

def debug_print_pretty(ds, desc, width=80):
    print "============== BEGIN {0} ================".format(desc)
    pprint.pprint(ds, width=width)
    print "============================================ END {0} ===".format(desc)

############
### Main ###
############

def main():

    ############################################################
    ### Stage 0a - Process the config file and start logging ###
    ############################################################

    # Note:
    try:
        _conf = Configuration()
        if 'DB_PROXY' in globals():
            _conf.configure(DB_PROXY)
        else:
            _conf.configure()
    except ConfigurationException as e:
        # Just re-raise normal ConfigurationExceptions
        raise e
    except Exception as e:
        # Somehow we got a Configuration section exception that was not a
        # ConfigurationException.  Since logging is not yet setup, we'll simply
        # print out the exception as normal and then exit.  This prevents
        # duplicate stacktraces from being printed.
        sys.stderr.write(traceback.format_exc())
        exit(1)

    # Set the global verbosity
    VERBOSE = _conf.general.verbose

    # Setup logging
    setup_logging(_conf, time.time())

    ###############################################
    ### Stage 0b - Read the sample fasta files. ###
    ###############################################

    # List of sample reads
    _samples = Samples()
    _samples._read_samples_filelist(_conf.general.samples_filelist)

    #################################################
    ### Stage 1 - Process the samples fasta files ###
    #################################################

    # Stores the sequences and their meta data
    _short_seqs = ShortSeqs(_conf)

    # Process
    with file_caching(_conf, 1, _short_seqs) as cached_objs:
        if not cached_objs:
            _short_seqs.process_sample_fastas(_samples)
        else:
            _short_seqs = cached_objs[0]

    #################################
    ### Stage 2 - Load from cache ###
    #################################
    if _conf.general.cache:
        _short_seqs.load_from_cache()

    ################################################
    ### Stage 3 - Generate the FASTA search file ###
    ################################################

    # First, confirm that _short_seqs is not empty.
    if not _short_seqs:
        raise CannotContinueException, \
            (   "Cannot continue: "
                "All sequence data have been filtered out by Prost.\n"
                "Try relaxing criteria in prost.config.\n"
                "For example, try lowering min_seq_count.")
    _short_seqs.generate_fasta_search_file(_conf.general.cache,
            _conf.general.input_file_for_alignments)

    ##########################################################
    ### Stage 4 - Perform Alignments (e.g. BBMap searches) ###
    ##########################################################
    _alignments = Alignments(_conf)
    _alignments.align(_conf)

    #########################################################
    ### Stage 5 - Designate Hits (Step One) and normalize ###
    #########################################################

    with file_caching(_conf, 5, _short_seqs) as cached_objs:
        if not cached_objs:
            _short_seqs.designation_step_one(_alignments.genome_alignment,
                    _conf.general.max_locations_allowed)
            _short_seqs.normalize()
        else:
            _short_seqs = cached_objs[0]

    ############################################
    ### Stage 6 - Caching - Writing to Cache ###
    ############################################
    if _conf.general.cache:
        _alignments.genome_alignment.write_to_cache(_short_seqs, DB_PROXY)

    ###########################################
    ### Stage 7 - Designate Hits - Step Two ###
    ###########################################
    _short_seqs.designation_step_two()

    #############################################
    ### Stage 8 - Binning by Genomic Location ###
    #############################################

    _gen_loc_bins = GenLocBins(_short_seqs)
    # Perform the binning
    _gen_loc_bins.perform_binning()
    # Sum counts across samples and save
    _gen_loc_bins.manipulate_per_sample_counts()
    _gen_loc_bins.post_binning_processing(_conf)

    ############################
    ### Stage 9 - Annotation ###
    ############################

    # I'm not sure why file caching ever worked with Annotation classes, because
    # thoses class have Enums as embedded classes, which Python 2.7 / Pickle
    # versions < 4 don't support.
    #
    # Because of this limitation, we may never support file caching from this
    # point onwards.
    #
    # Disabling file caching here for now.

    progress = Progress("Annotation")
    logging.info("Annotation...")
    sys.stderr.write('\n')
    _short_seqs.perform_annotations(_alignments, _conf.general.species)
    progress.done()

    #########################################
    #### Stage 10 - Binning by Annotation ###
    #########################################

    _annotation_bins = AnnotationBins(_short_seqs, _gen_loc_bins)
    progress = Progress("Binning by Annotation")
    _annotation_bins.perform_binning()
    # Sum counts across samples and save
    _annotation_bins.manipulate_per_sample_counts()
    _annotation_bins.post_binning_processing()
    progress.done()

    ##################################
    ### Stage 11 - Binning by Seed ###
    ##################################

    _seed_bins = SeedBins(_short_seqs)
    _seed_bins.perform_binning(_annotation_bins)
    # Sum counts across samples and save
    _seed_bins.manipulate_per_sample_counts()

    #############################
    ### Stage 12 - Output TSV ###
    #############################

    Output().tsv_isomirs(_conf, _samples, _short_seqs, _alignments,
            _gen_loc_bins, _seed_bins, _annotation_bins)
    Output().tsv_comp_by_gen_loc(_conf, _samples, _short_seqs, _gen_loc_bins,
            _alignments)
    Output().tsv_comp_by_annotation(_conf, _samples, _short_seqs,
            _annotation_bins, _alignments)
    Output().tsv_comp_by_seed(_conf, _samples, _short_seqs, _seed_bins,
            _gen_loc_bins, _alignments)
    Output().tsv_mirror_mirs(_conf, _short_seqs, _gen_loc_bins, _alignments)
    Output().tsv_arm_switch(_conf, _samples, _short_seqs, _annotation_bins, _alignments)
    Output().tsv_no_hits(_conf, _samples, _short_seqs, _alignments)

    ###############################
    ### Stage 13 - Output Excel ###
    ###############################

    # TODO: Provide user-configurable options to user to allow user to turn
    # on/off --one mode and --many mode.  For now, just go.
    progress = Progress("Generating Excel Spreadsheet")
    prefix = _conf.general.output_file_prefix
    prost.excel.create_one(prefix, prefix, _conf.general.output_files)
    progress.done()


# vim: softtabstop=4:shiftwidth=4:expandtab
