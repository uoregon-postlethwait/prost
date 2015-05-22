.. _index:

********
*Prost!*
********

*Prost!* (PRocessing Of Small Transcripts) is a python application that
quantifies and annotates microRNA (miRNA) expression in chordates and
vertebrates with assembled genomes.  *Prost!* works by counting short
transcripts within a user-specifiable length range.  These counted transcripts
are aligned to a user specifiable genome allowing for post-transcriptional
modification (e.g. untemplated additions, editing, alternative cutting) and
then "binned" together based on genomic location.  Each bin is then annotated
with miRBase mature sequences and hairpins, as well as other types of RNA
obtained from Ensembl's Biomart.

Site Contents
=============

In addition to this home page, there are several other pages on this site as
well. Click on the *Site* link in the header to access those pages.  Here is an
overview of the content of the other pages.

.. toctree::
   :maxdepth: 2

   demo
   output
   bbmap
   biomart

..
  Add a hidden table of contents that includes the development page, so
  that the 'Site' link in the header includes the Development page.  We don't
  wish to include it in the Site Contents as well.

.. toctree::
   :hidden: 
   
   development

.. _prost_installation:

Installation of *Prost!*
========================

Installation via pip is probably the easiest way to get *Prost!* installed::

    pip install prost

Alternatively, you can download a release and install directly with setup.py::
    
    # Point your browser at 
    https://github.com/uoregon-postlethwait/prost/releases/latest 

    # And use your browser to download the latest release.  For example:
    https://github.com/uoregon-postlethwait/prost/archive/v0.7.3.tar.gz
    
    # Extract
    tar xzf prost-0.7.3.tar.gz
    cd prost-0.7.3
    
    # Install (systemwide)
    python setup.py install
    
    # Install (in your home directory)
    python setup.py install --user

Finally, you can clone the git repository and install directly with setup.py::

    # Clone repository
    git clone https://github.com/uoregon-postlethwait/prost
    cd prost

    # Install (systemwide)
    python setup.py install
    
    # Install (in your home directory)
    python setup.py install --user

Quickstart
==========

To help you get *Prost!* up and running quickly, we've provided a demo in the
style of a tutorial.  The demo includes all the sample data you'll need (except
for the reference genome due to size constraints, which you'll have to download
separately).

To get started, see the :ref:`demo` page.

Requirements and Dependencies
=============================

*Prost!* is designed to be easily run from the command line.  A basic familiarity
of the command line and optionally the ability to edit a configuration file are
all that is required to *Prost!*.

Listed below are the various requirements that need to be fulfilled in order to
run *Prost!*.

Hello

.. _system_and_software_dependencies: 

System and software dependencies
````````````````````````````````

* A Linux or Mac OS X environment.
* Python version 2.7.x (tested with Python 2.7.9, 2.7.5, and 2.7.2).
* BBMap short read aligner, available for download at 
  `SourceForge <http://sourceforge.net/projects/bbmap/>`_.  See
  :ref:`bbmap_installation`.

Configuration files
```````````````````

* The *Prost!* configuration file (default ``prost.config``). This is
  a basic configuration file that is well documented and should be easy to
  quickly comprehend and edit.
* A list of high throughput sequencing FASTA format files of samples
  (default ``samples_filelist``). The format is ``fileName descriptiveName``,
  with one samples file per line.  

Samples FASTA files
```````````````````

* Preprocessed input fasta files containing short sequence reads from your
  samples.  Quality filtered and barcodes or adapters removed so that the only
  sequences remaining in the files are full, high quality, short sequences.

BBMap databases
```````````````

Four BBMap databases are required by *Prost!*.  See :ref:`building_bbmap_dbs`
for instructions on building a BBMap database.

The following databases are required by *Prost*!:

* A reference genome BBMap database.  The reference genome does not need to be
  annotated or assembled to high N50s.  We have tested *Prost!* on PE250
  Illumina sequencing reads on a non-model organism with good results.
* Annotation BBMap databases.  NOTE: At this time, nucleotides uracil and
  thymine must be coded as **T**'s (and *not* as **U**'s) in the annotation
  FASTA files.

  * A database of annotated mature microRNAs. This is usually built directly from
    `miRBase's "mature.fa" <ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz>`_ 
    (or from an augmented version of that FASTA file), and should contain
    mature miRNAs from several species (including the species you are
    studying).  The FASTA header file must follow miRBase's convention of
    prefixing each miRNA with a three letter species abbreviation and a dash.
    For example::

        >dre-miR-451 MIMAT0001634 Danio rerio miR-451
        AAACCGTTACCATTACTGAGTT
  * A database of annotated microRNA hairpins. This is usually built directly from
    `miRBase's "hairpin.fa" <ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz>`_ 
    (or from an augmented version of that FASTA file), following the same
    conventions for the mature miRNAs described above.
  * A database of annotated "other RNAs" (e.g snoRNAs and lincRNAs) for the
    species under study. See :doc:`biomart`.  Following those instructions will
    produce correctly formatted FASTA files (the FASTA header format is:
    ``>geneName|biotype|geneID``).


*Prost!* Output
```````````````

See :doc:`output` for a description of the output of *Prost!*.

.. _command_line_help:

Command Line Help
=================

*Prost!* provides configuration options provided through the command line.  To
see a list of available options with contextual help, as well as each options
defaults, type the following in a terminal window:

.. code-block:: bash

      prost -h
      
      # or
      prost --help

Configuration
=============

*Prost!* follows the usual convention of allowing configuration via defaults,
a configuration file, and command line flags.  Sensible defaults are provided,
so in general very few configuration options will be needed by the user.

A minimal configuration file (by default named ``prost.config``) is required by
*Prost!*.  At minumum, the configuration file needs to specify the species
under study (e.g. 'mmu' for mouse) in the **General** section.  In addition,
the **GenomeAlignment** section needs to specify the following:

* **name**: A user-defined name for the alignment (default: ``genome``).
* **tool**: The alignment tool being used (default: ``bbmap``; currently only
  BBMap is supported).
* **db**: The sequence database (e.g. a BBMap database) to be used for the Genome
  Alignment.
* **max_3p_mismatches**: The maximum number of mismatches allowed for each sequence
  alignment hit on the 3´ end of the hit (default: ``3``)
* **max_non_3p_mismatches**: The maximum number of mismatches allowed for each sequence
  alignment hit anywhere but the 3´ end of the hit (default: ``2``) (i.e.
  either on the 5´ end or in the middle of the sequence).
* **allow_indels**: Whether or not to allow alignments which have insertions or deletions.
* **indelnt_penalty_multiplier**: Multiply the number of inserted or deleted
  nucleotides in the alignment by this penalty, and add that to the number of
  3´ mismatches when determining whether to keep or reject that alignment.
  Basically, this can nearly eliminate indels from your dataset if that is your
  desire.

Below is an example of a minimal configuation file:

.. code-block:: ini

  [General]
  species: mmu
  samples_filelist: data/samples_filelist

  [GenomeAlignment]
  name: genome
  tool: bbmap
  db: /path/to/databases/Danio_rerio.Zv9.dna.toplevel
  max_3p_mismatches: 3
  max_non_3p_mismatches: 2
  allow_indels: yes

The configuration file as described above will not perform any annotations.
See the file ``prost.config`` for a working example of annotation alignment
sections.  The annotation alignment sections follow the same structure as the
**GenomeAlignment** section.

Currently, (nearly) every option in the **General** section of the
configuration file can also be controlled via command line flags.  See
:ref:`command_line_help`.

Running *Prost!* On Your Own Dataset
====================================

To run *Prost!* on your own dataset, you'll need to edit the configuration file.
In particular, you'll need to edit the fields 'species' and 'samples_filelist'
(in the **General** section) and the 'db' fields (in the **Alignment** sections).
Below is a snippet of the configuration file that shows you roughly what will
need to be edited (some additional fields are shown below for context, but do not
need to be edited):

.. code-block:: ini

   [General]
   species: dre
   samples_filelist: samples_filelist

   [GenomeAlignment]
   db: BBMap/Danio_rerio.Zv9.dna.toplevel

   [AnnotationAlignment1]
   type: MirbaseMirAnnotation
   db: BBMap/mature_miRBase21

   [AnnotationAlignment2]
   type: MirbaseHairpinAnnotation
   db: BBMap/hairpin_miRBase21

   [AnnotationAlignment3]
   type: BiomartOtherRNAAnnotation
   db: BBMap/BioMart_Dre79_otherRNA

After you have made those changes, simply run *Prost!* again:

.. code-block:: bash

   python prost

Funding
=======

*Prost!* has been funded by the following grants:

* Identification of MiRNAs Involving Midfacial Development and Clefting; NIH - National Institute of Dental and Craniofacial Research (U01 DE020076)
* Advancing the Scientific Potential of Transcriptomics in Aquatic Models; NIH - Office of the Director (R24 OD011199)
* Resources for Teleost Gene Duplicates and Human Disease; NIH - Office of the Director (R01 OD011116)
* Mechanisms of Sex Determination in Zebrafish; NIH - National Institute of General Medical Sciences (R01 GM085318)
* Developmental Mechanisms for the Evolution of Bone Loss; NIH - National Institute on Aging (R01 AG031922)
* Signaling Hierarchies in Vertebrate Development: CP1:  A zebrafish model of phenotypic variation associated with Fraser syndrome; NIH - Eunice Kennedy Shriver National Institute of Child Health and Human Development (P01 HD22486)

.. Hyperlinks
.. _Python: http://www.python.org/
.. _mature.fa: ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
.. _hairpin.fa: ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz

.. If you want 
   .. automodule:: prost
   :members:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

