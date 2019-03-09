.. _index:

********
*Prost!*
********

*Prost!* (PRocessing Of Small Transcripts) is a python application that
quantifies and annotates microRNA (miRNA) expression in metazoans with
assembled genomes.  *Prost!* works by counting short transcripts within a
user-specifiable length range.  These counted transcripts are aligned to a user
specifiable genome allowing for post-transcriptional modification (e.g.
untemplated additions, editing, alternative cutting) and then "binned" together
based on genomic location.  Each bin is then annotated with databases of mature
miRNAs, hairpins, and other types of RNAs (the databases may be derived from
miRBase, Ensembl's BioMart, other databases, or may be custom built by the
user).

Authors
=======

* Peter Batzel - Conceptual Design, Software Engineering, and Algorithm Design
* Thomas Desvignes - Conceptual Design Leader
* Jason Sydes - Conceptual Design, Software Engineering, and Algorithm Design
* Brian F. Eames - Prototype Design
* John H. Postlethwait - Project Advisor

News
====

**March 8th, 2019**: *Prost!* was published in 
`Scientific Reports <https://www.nature.com/articles/s41598-019-40361-8>`_.

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

Installation via pip is probably the easiest way to get *Prost!* installed. 
To install to your home directory::

    pip install prost --user

Or you can install systemwide (if you have sudo privileges)::

    sudo pip install prost

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

* The *Prost!* configuration file (default ``prost.config``).  We recommend
  you download our example *Prost!* configuration file 
  (`prost.config.example <https://raw.githubusercontent.com/uoregon-postlethwait/prost/master/prost.config.example>`_)
  and modify it to suit your experiments.  The example configuration file is
  well documented and should be easy to quickly comprehend and edit.
* A list of high throughput sequencing FASTA format files of samples
  (default ``samples_filelist``). The format is ``fileName descriptiveName``,
  with one samples file per line.  

Samples FASTA files
```````````````````

* Preprocessed input fasta files containing short sequence reads from your
  samples.  Quality filtered and barcodes or adapters removed so that the only
  sequences remaining in the files are full, high quality, short sequences.

Reference genome BBMap database
```````````````````````````````

*Prost!* requires a BBMap database of a reference genome of your species of
interest (or a closely related species). The reference genome does not need to
be annotated or assembled to high N50s. We have tested *Prost!* on PE250
Illumina sequencing reads on a non-model organism with good results. See
:ref:`building_bbmap_dbs` for instructions on building a BBMap database.

Annotation FASTA files
``````````````````````

Three annotation FASTA files are required by *Prost!*. (NOTE: At this time,
nucleotides uracil and thymine must be coded as **T**'s (and *not* as **U**'s)
in the annotation FASTA files.)

* A FASTA file of annotated mature microRNAs. This is usually built directly from
  `miRBase's "mature.fa" <ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz>`_ 
  (or from an augmented version of that FASTA file), and should contain
  mature miRNAs from several species (including the species you are
  studying).  The FASTA header file must follow miRBase's convention of
  prefixing each miRNA with a three letter species abbreviation and a dash.
  For example::

      >dre-miR-451 MIMAT0001634 Danio rerio miR-451
      AAACCGTTACCATTACTGAGTT
* A FASTA file of annotated microRNA hairpins. This is usually built directly from
  `miRBase's "hairpin.fa" <ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz>`_ 
  (or from an augmented version of that FASTA file), following the same
  conventions for the mature miRNAs described above.
* A FASTA file of annotated "other RNAs" (e.g snoRNAs and lincRNAs) for the
  species under study, or a related species.   See :doc:`biomart`.  Following
  those instructions will produce correctly formatted FASTA files (the FASTA
  header format is: ``>geneName|biotype|geneID``). *Prost!* requires a BioMart
  file because it usually provides useful annotations.  If you do not wish to
  supply a real BioMart file, simply download and point *Prost!* at the
  following fake BioMart file: `fake_biomart_file.fa <https://raw.githubusercontent.com/uoregon-postlethwait/prost/master/fake_biomart_file.fa>`_

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
*Prost!* (you can download and modify our example configuration file here: 
`prost.config.example <https://raw.githubusercontent.com/uoregon-postlethwait/prost/master/prost.config.example>`_). 
At minumum, the configuration file needs to specify the species
under study (e.g. 'mmu' for mouse) in the **General** section.  In addition,
the **GenomeAlignment** section needs to specify the following:

* **name**: A user-defined name for the alignment (default: ``genome``).
* **tool**: The alignment tool being used (default: ``bbmap``; currently only
  BBMap is supported).
* **db**: The BBMap database of the reference genome (e.g. a BBMap database) to
  be used for the Genome Alignment.  (Note that for the AnnotationAlignment sections below,
  this **db** field must point to a FASTA file instead of a BBMap database!)
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
  db: /path/to/databases/Danio_rerio.GRCz10.dna.toplevel
  max_3p_mismatches: 3
  max_non_3p_mismatches: 2
  allow_indels: yes

The configuration file as described above will not perform any annotations.
See the file `prost.config.example <https://raw.githubusercontent.com/uoregon-postlethwait/prost/master/prost.config.example>`_ 
for a working example of annotation alignment sections.  The
**AnnotationAlignment** sections follow the same structure as the
**GenomeAlignment** section; the only difference being that the **db** field
must point to a FASTA file instead.

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
   mature_mir_annotation_fasta: BBMap/mature_miRBase21.fa

   [GenomeAlignment]
   db: BBMap/Danio_rerio.GRCz10.dna.toplevel

   [AnnotationAlignment1]
   type: MirbaseMirAnnotation
   db: BBMap/mature_miRBase21.fa

   [AnnotationAlignment2]
   type: MirbaseHairpinAnnotation
   db: BBMap/hairpin_miRBase21.fa

   [AnnotationAlignment3]
   type: BiomartOtherRNAAnnotation
   db: BBMap/BioMart_Dre89_otherRNA.fa

After you have made those changes, simply run *Prost!* again:

.. code-block:: bash

   prost

Funding
=======

*Prost!* has been funded by the following grants:

* Identification of MiRNAs Involving Midfacial Development and Clefting; NIH - National Institute of Dental and Craniofacial Research (U01 DE020076)
* Advancing the Scientific Potential of Transcriptomics in Aquatic Models; NIH - Office of the Director (R24 OD011199)
* Resources for Teleost Gene Duplicates and Human Disease; NIH - Office of the Director (R01 OD011116)
* Mechanisms of Sex Determination in Zebrafish; NIH - National Institute of General Medical Sciences (R01 GM085318)
* Developmental Mechanisms for the Evolution of Bone Loss; NIH - National Institute on Aging (R01 AG031922)
* Signaling Hierarchies in Vertebrate Development: CP1:  A zebrafish model of phenotypic variation associated with Fraser syndrome; NIH - Eunice Kennedy Shriver National Institute of Child Health and Human Development (P01 HD22486)
* Antarctic Fish and MicroRNA Control of Development and Physiology; NSF - Office of Polar Program (OPP #154338)

Citing *Prost!*
===============

Thomas Desvignes, Peter Batzel, Jason Sydes, B. Frank Eames, and John H. Postlethwait. 2019. “MiRNA Analysis with Prost! Reveals Evolutionary Conservation of Organ-Enriched Expression and Post-Transcriptional Modifications in Three-Spined Stickleback and Zebrafish.” Scientific Reports 9 (1): 3913. https://doi.org/10.1038/s41598-019-40361-8.

.. Hyperlinks
.. _Python: http://www.python.org/
.. _mature.fa: ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
.. _hairpin.fa: ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz

.. If you want 
   .. automodule:: prost
   :members:
