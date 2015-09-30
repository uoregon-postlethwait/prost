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

################
### Defaults ###
################

# Whether or not to be verbose
VERBOSE =                                                              True

# The configuration file
DEFAULT_PROST_CONFIG_FILE =                                  "prost.config"

# The default output directory
DEFAULT_PROST_OUTPUT_DIR =                                   "prost.output"

# The default samples filelist
DEFAULT_SAMPLES_FILELIST =                               "samples_filelist"

# The default search input file for the genome and annotation alignments
DEFAULT_SEARCH_INPUT_FILE_FOR_ALIGNMENTS =                      "search.fa"

# The default output file prefix
DEFAULT_OUTPUT_FILE_PREFIX =                                  "prost_output"

# The default minimun number of times a sequence must appear (in total) across
# all samples to be included in the data set.
DEFAULT_MIN_SEQ_COUNT =                                                  30

# The amount of genomic wiggle
DEFAULT_WIGGLE =                                                          5

# Minimum and maximum sequence length allowed for a short sequence.
DEFAULT_MIN_SEQ_LENGTH =                                                 14
DEFAULT_MAX_SEQ_LENGTH =                                                 28

# The maximum number of threads (processors) to use for various operations.
DEFAULT_MAX_THREADS =                                                     4

# Default maximum number of locations to report
DEFAULT_MAX_LOCATIONS_TO_REPORT =                                        20

# Default maximum number of locations allowed per sequence
DEFAULT_MAX_LOCATIONS_ALLOWED =                                          40

# The size of buckets (in nucleotides) for mirror-miR detection, and the amount
# of overlap (in nucleotides) for a mirror-miR.  Note that if you're using
# negative overlap, you'll want the bucket length to be some X times larger than
# the absolute value of the overlap length.  Haven't done the math, but 6 times
# larger should do. If doing positive overlap, then 200bp for the bucket length
# should be fine.
DEFAULT_MIRROR_BUCKET_LENGTH =                                          200
DEFAULT_MIRROR_OVERLAP =                                                  5
#DEFAULT_MIRROR_BUCKET_LENGTH =                                         1200
#DEFAULT_MIRROR_OVERLAP =                                               -200
#DEFAULT_MIRROR_BUCKET_LENGTH =                                        12000
#DEFAULT_MIRROR_OVERLAP =                                              -2000

# Default value for alignment configuration option
# 'indelnt_penalty_multiplier'.
DEFAULT_INDELNT_PENALTY_MULTIPLIER =                                      3


#################
### Constants ###
#################

# File suffixes for output files.
OUTPUT_FILE_SUFFIX_ISOMIRS =                         "uncompressed_isomiRs"
OUTPUT_FILE_SUFFIX_NO_HITS =                              "no_genomic_hits"
OUTPUT_FILE_SUFFIX_COMP_BY_GEN_LOC =       "compressed_by_genomic_location"
OUTPUT_FILE_SUFFIX_COMP_BY_SEED =                      "compressed_by_seed"
OUTPUT_FILE_SUFFIX_COMP_BY_ANNOTATION =          "compressed_by_annotation"
OUTPUT_FILE_SUFFIX_MIRROR_MIRS =                    "candidate_mirror-mirs"
OUTPUT_FILE_SUFFIX_ARM_SWITCH =                      "candidate_arm_switch"

# Collect all output file suffixes for easier machine parsing
OUTPUT_FILE_SUFFIXES = [            OUTPUT_FILE_SUFFIX_ISOMIRS,
                                    OUTPUT_FILE_SUFFIX_COMP_BY_GEN_LOC,
                                    OUTPUT_FILE_SUFFIX_COMP_BY_SEED,
                                    OUTPUT_FILE_SUFFIX_COMP_BY_ANNOTATION,
                                    OUTPUT_FILE_SUFFIX_MIRROR_MIRS,
                                    OUTPUT_FILE_SUFFIX_ARM_SWITCH,
                                    OUTPUT_FILE_SUFFIX_NO_HITS]

# The max length of the 'designtation' varchar column in the genomic_location
# db table.
COL_LEN__GENOMIC_LOCATION__DESIGNATION =                                  6

# The designation for a 'one' sequence.
DESIGNATION_ONE =                                                         1

# The designation for a 'two' sequence.
DESIGNATION_TWO =                                                         2

# The designation for a 'two' sequence.
DESIGNATION_THREE =                                                       3

# A magic number representing a sequence designation of either 'one' or 'two'.
DESIGNATION_ONE_OR_TWO =                                                -12

# A magic number representing a sequence designation for sequences that have
# exceeded max_locations_allowed.
DESIGNATION_MAX_LOCATIONS_ALLOWED_EXCEEDED =                             98

# A magic number representing a sequence designation of 'no hits in genome db'.
DESIGNATION_NO_HIT =                                                     99

# TODO: rename - This has to do with illumina reads
ONE_MILLION =                                                     1000000.0

# Annotation Indexes
SPECIES_IDX =                                                             0
OTHER_SPECIES_IDX =                                                       1

# For now, limit the number of genomic locations that a bin can have if it
# is to be considered for mirror mir detection.
TODO__HARDCODE_SMALL_MAX_GEN_LOC_LEN_FOR_MIRROR_MIR_DETECTION = 3

# Hardcoded BBMap minid parameter.
TODO_HARDCODED_MINID = "0.50"

# This is the suffix added to the designation to indicate that it has an
# insertion or deletion.
INDEL_INDICATOR =                                                       "g"

# Executables
EXECUTABLE_PYTHON =                                                'python'
EXECUTABLE_BBMAP =                                        'bbmapskimmer.sh'

# Prost log filename
LOG_FILENAME =                                                  'prost.log'
# Overall log level
LOG_LEVEL =                                                         'DEBUG'
# Log level for prost.log
LOG_LEVEL_FOR_FILE =                                                 'INFO'
# Logging configuration dict
LOGGING = {
    'version': 1,
    'handlers': {
        #'console': {
        #    'class': 'logging.StreamHandler',
        #    'level': LOG_LEVEL,
        #},
        'file': {
            'class': 'logging.FileHandler',
            'level': LOG_LEVEL_FOR_FILE,
            'filename': LOG_FILENAME,
            'mode': 'w',
        }
    },
    'root': {
        'level': LOG_LEVEL,
        #'handlers': ['console', 'file']
        'handlers': ['file']
    },
}

# Note:
# * CONFIG_FILE* : These constants refer to the config file only.
# * CONFIGURATION* : These constants refer to the combination of the settings
#       found in the config file and command line arguments.

# Config file specifications
CONFIG_FILE_REQUIRED_GENERAL_FIELDS =                           ['species']
CONFIG_FILE_REQUIRED_ALIGNMENT_FIELDS = \
                                               ['name',
                                                'tool',
                                                'db',
                                                'max_3p_mismatches',
                                                'max_non_3p_mismatches',
                                                'allow_indels']

CONFIG_FILE_ADDITIONAL_REQUIRED_ANNOTATION_ALIGNMENT_FIELDS = \
                                                                   ['type']
CONFIG_FILE_SUPPORTED_ALIGNMENT_TOOLS =                           ['bbmap']
CONFIG_FILE_SUPPORTED_ANNOTATION_ALIGNMENT_TYPES = ['MirbaseMirAnnotation',
                                                'MirbaseHairpinAnnotation',
                                               'BiomartOtherRNAAnnotation']
CONFIG_FILE_SUPPORTED_DATABASES =                                 ['mysql']
CONFIG_FILE_REQUIRED_CACHE_FIELDS = \
                               ['type', 'database', 'username', 'password']
CONFIG_FILE_BOOLEAN_FIELDS =                   ['verbose',
                                                'read_from_file_cache',
                                                'write_to_file_cache',
                                                'skip_sequence_alignments',
                                                'skip_genome_alignment',
                                                'skip_annotation_alignments',
                                                'cache',
                                                'cache_read',
                                                'cache_write',
                                                'create_tables',
                                                'allow_indels']
CONFIG_FILE_INT_FIELDS =                       ['max_3p_mismatches',
                                                'max_non_3p_mismatches',
                                                'indelnt_penalty_multiplier',
                                                'min_seq_count',
                                                'min_seq_length',
                                                'max_seq_length',
                                                'max_threads',
                                                'max_locations_to_report',
                                                'max_locations_allowed',
                                                'wiggle']
# Configuration fields that should be non-negative
CONFIGURATION_NON_NEGATIVE_FIELDS =            ['max_3p_mismatches',
                                                'max_non_3p_mismatches',
                                                'wiggle']

# Configuration fields that should be positive
CONFIGURATION_POSITIVE_FIELDS =                ['indelnt_penalty_multiplier',
                                                'min_seq_count',
                                                'min_seq_length',
                                                'max_seq_length',
                                                'max_threads',
                                                'max_locations_to_report',
                                                'max_locations_allowed']


# BBMap default parameters
# Generate MD tags
BBMAP_OPT_DEFAULT_MDTAG =                                               "t"
# Generate cigar strings (defaults to sam 1.4 extended cigar format)
BBMAP_OPT_DEFAULT_CIGAR =                                               "t"
# Helpful tags if debugging:
BBMAP_OPT_DEFAULT_SCORETAG =                                            "f"
BBMAP_OPT_DEFAULT_STOPTAG =                                             "f"
BBMAP_OPT_DEFAULT_IDTAG =                                               "f"
BBMAP_OPT_DEFAULT_INSERTTAG =                                           "f"
# Do alignments faster, in C code. Requires compiling the C code; details are
# in /jni/README.txt.
# TODO: Make this a prost configurable option.
BBMAP_OPT_DEFAULT_USEJNI =                                              "f"
# Set to false if unmapped reads should not be printed to 'out=' target (saves
# time and disk space).
BBMAP_OPT_DEFAULT_OUTPUTUNMAPPED =                                      'f'
# Truncate read and ref names at the first whitespace, assuming that the
# remainder is a comment or description.
BBMAP_OPT_DEFAULT_TRIMREADDESCRIPTIONS =                                "t"
# Kmer length of index.  Lower is slower, more sensitive, and uses slightly
# less RAM.  Range is 7 to 15.
BBMAP_OPT_DEFAULT_K =                                                   "7"
# Throw away ~80% of kmers based on remainder modulo a number (reduces RAM by
# 50% and sensitivity slightly). Should be enabled both when building the index
# AND when mapping.
BBMAP_OPT_DEFAULT_USEMODULO =                                           "f"
# When enabled, do not allow indels longer than 'maxindel'. By default these
# are not sought, but may be found anyway.
BBMAP_OPT_DEFAULT_STRICTMAXINDEL =                                      "t"
# Set behavior on ambiguously-mapped reads (with multiple top-scoring mapping
# locations) to 'all':
#       multiple top-scoring mapping locations).
#       best    (use the first best site)
#       toss    (consider unmapped)
#       random  (select one top-scoring site randomly)
#       all     (retain all top-scoring sites)
BBMAP_OPT_DEFAULT_AMBIGUOUS =                                         "all"
# Do not print secondary alignments
BBMAP_OPT_DEFAULT_SECONDARY =                                           "t"
# Print only secondary alignments with score of at least this fraction of
# primary. Default 0.95.
BBMAP_OPT_DEFAULT_SSSR =                                             "0.25"
# (secondarysiteasambiguousonly) Only print secondary alignments for
# ambiguously-mapped reads.
BBMAP_OPT_DEFAULT_SSAO =                                                "f"
# Maximum number of total alignments to print per read. Only relevant when
# secondary=t.
BBMAP_OPT_DEFAULT_MAXSITES =                                      "4000000"
# Don't analyze (or print) more than this many alignments per read.
BBMAP_OPT_DEFAULT_MAXSITES2 =                                     "4000000"
# This flag is a macro which sets other paramters to run slower, at greater
# sensitivity.  'vslow' is even slower.
BBMAP_OPT_DEFAULT_SLOW =                                                "t"

# Turn on additional debug timing.
DEBUG_TIMING =                                                        False

# vim: softtabstop=4:shiftwidth=4:expandtab
