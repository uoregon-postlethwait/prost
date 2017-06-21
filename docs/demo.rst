.. _demo:

****************
Quick Start Demo
****************

This demo and tutorial will walk you through getting *Prost!* running on a
provided zebrafish test dataset.

Requirements
============

*Prost!* is designed to be easily run from the command line.  A basic familiarity
of the command line and optionally the ability to edit a configuration file are
all that is required to *Prost!*.

Before continuing, make sure you've met the :ref:`system_and_software_dependencies`.

Dependencies Not Provided in Demo
---------------------------------

Due to size constraints, a BBMap database of the zebrafish reference genome is
not included.  Instructions on how to download and build a BBMap database of
the zebrafish genome are included below.

Dependencies Provided in Demo
-----------------------------

This demo provides several dependencies via the GitHub repository
https://github.com/uoregon-postlethwait/prost-demo. If you wish to run
*Prost!* against a different dataset, you'll likely need to provide some or all
of these requirements.  The following dependencies are provided:

* An appropriate *Prost!* configuration file (``prost.config``).
* A list of high throughput sequencing fasta format files of samples
  (``samples_filelist``) and the corresponding samples files
  (``samples/sample.brain1.fa``, ``samples/sample.heart2.fa``, and
  ``samples/sample.ovary3.fa``).
* Annotation FASTA files:

  * A FASTA file of mature microRNAs in miRBase release 21
    (``fa/mature_miRBase21.fa.gz``).
  * A FASTA file of microRNA hairpins in miRBase release 21
    (``fa/hairpin_miRBase21.fa.gz``).
  * A FASTA file of "other RNAs" (e.g. snoRNAs) found in zebrafish, provided by
    Ensembl BioMart release 79 (``fa/BioMart_Dre79_otherRNA.fa.gz``).
* A script (``scripts/setup.sh``) which uncompresses the sample and annotation
  FASTA files described above.

Quick Start
===========

#. Make sure you've met the following requirements:

   * :ref:`system_and_software_dependencies`.
   * You have *Prost!* installed (see :ref:`prost_installation`)
   * You have `Java 7 <http://www.oracle.com/technetwork/java/javase/downloads/jre7-downloads-1880261.html>`_ installed (and enabled).
   * You have BBMap installed (see :ref:`bbmap_installation`)
#. If you have a slow internet connection, you may wish to begin downloading
   the zebrafish reference genome (see step 5).
#. Download the *Prost!* demo. 

   * You can do this manually by navigating your browser to the prost-demo
     Releases page and clicking on the `latest release link <https://github.com/uoregon-postlethwait/prost-demo/releases/latest>`_.
   * Or you do it automatically with the following command::

           curl -SLOJ $(curl -s https://api.github.com/repos/uoregon-postlethwait/prost-demo/releases/latest | grep zipball_url | head -n 1 | cut -d\" -f4)
#. Extract the *Prost!* demo and cd into the demo directory::

        # Replace X.X.X with the version of the demo.
        unzip prost-demo-X.X.X.zip
        cd prost-demo-X.X.X
#. Run the prost-demo's setup script.  This script uncompresses the samples
   files as well as builds the annotation BBMap databases.  Run the setup 
   script like so::

        sh scripts/setup.sh
#. Download and extract the zebrafish reference genome.  (Alternatively, if you
   already have a copy of the zebrafish assembly available, you may use that
   instead.  It doesn't need to be the exact version shown below.  Just give
   ``bbmap.sh`` the full PATH to the ``ref=`` argument in the following step.)::

        # Using wget:
        wget ftp://ftp.ensembl.org/pub/release-89/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.toplevel.fa.gz

        # Or using curl:
        curl -O ftp://ftp.ensembl.org/pub/release-89/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.toplevel.fa.gz
        
#. Build a BBMap database of the zebrafish reference genome::

        # We've found a k-mer length of 7 or 8 works best for these data:
        bbmap.sh k=7 path=BBMap/Danio_rerio.GRCz10.dna.toplevel ref=Danio_rerio.GRCz10.dna.toplevel.fa.gz

#. Run *Prost!*::

        # The demo defaults to using 4 threads.
        prost

        # Alternatively, you may run Prost! on more CPUs if they are available.
        prost --max-threads 12

Prost Logging Output in Demo
============================

If everything went well, *Prost!* will output the following to both the console
(i.e. STDOUT), as well as to ``prost.log``, and it should look something like
this::

        ∴ prost --max-threads 12
        Prost! version 0.7.3.                            5:00:41 PM PDT on Apr 16, 2015.

        Reading in all fasta files...
        reading file1 (sample.brain1.fa)...
            83157/83157 so far... done. (elapsed time: 0.2s)
        reading file2 (sample.heart2.fa)...
            24039/24039 so far... done. (elapsed time: 0.1s)
        reading file3 (sample.ovary3.fa)...
            18492/18492 so far... done. (elapsed time: 0.1s)
        Rejecting low reads seqs... done. (elapsed time: 0.0s)
        Writing fasta search file... done. (elapsed time: 0.0s)
        Alignments... done. (elapsed time: 77.9s (or 1m17s))
        Designation step ONE...
                Reading BBMap hits from file...
                    244/244 so far... done. (elapsed time: 0.0s)
        Designation step ONE... done. (elapsed time: 0.0s)
        Normalization: calculating per-sample totals.... done. (elapsed time: 0.0s)
        Normalization: normalizing read counts.... done. (elapsed time: 0.0s)
        Designation step TWO...
            14/14 so far... done. (elapsed time: 0.0s)
        Annotation...
                Reading BBMap hits from file...
                    1113/1113 so far... done. (elapsed time: 0.1s)
                Reading BBMap hits from file...
                    2478/2478 so far... done. (elapsed time: 0.1s)
                Reading BBMap hits from file...
                    238/238 so far... done. (elapsed time: 0.0s)
        Annotation... done. (elapsed time: 0.2s)
        Binning by Annotation... done. (elapsed time: 0.0s)
        Writing output file...
            65/65 so far... done. (elapsed time: 0.1s)
        Writing comressed output file...
            4/4 so far... done. (elapsed time: 0.1s)
        Writing annotation-comressed output file...
            0 so far... done. (elapsed time: 0.1s)
        Writing seed-comressed output file...
            0 so far... done. (elapsed time: 0.1s)
        Writing mirror miRs output file... done. (elapsed time: 0.1s)
        Writing arm switch output file... done. (elapsed time: 0.1s)
        Writing no_hits output file...
            65/65 so far... done. (elapsed time: 0.1s)
        Generating Excel Spreadsheet... done. (elapsed time: 0.2s)
        Total Prost running time: 81.2s (or 1m21s).


Prost Excel Output
==================

Once complete, *Prost!* produces an Excel spreadsheet with several tabs (as
well as several tab separated value (TSV) files which are identical to the
Excel tabs, minus formatting).  The by_gen_loc bin is a good place to start.
Please see the documentation for descriptions of each tab and column.  

.. todo: (above) Add in a link to that documentation when it exists.

.. todo: (below) Not sure, this section is duplicated...

Adapting the Demo for Your Dataset
==================================

You can adapt this demo to easily run *Prost!* on your own dataset.  To do so,
you'll need to edit the configuration file ``prost.config``. In particular,
you'll need to edit these fields in the **General** section:

* *species* - to specify your species
* *samples_filelist* - to point to your file list of samples;
  alternatively, you can simply edit the file ``samples_filelist`` in the
  current directory.

You may also need to edit the *db* fields in the **Alignment** sections if you
are using different genome or annotation databases.

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
   db: BBMap/BioMart_Dre79_otherRNA.fa

After you have made those changes, simply run *Prost!* again:

.. code-block:: bash

   prost

.. If you enable this, it will break your left TOC in the RTD theme.
   .. toctree::
   :maxdepth: 2
.. Instead, ... you might have to manually specify the TOC...

.. Hyperlinks
.. _Python: http://www.python.org/
.. _mature.fa: ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
.. _hairpin.fa: ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz

.. If you want 
   .. automodule:: prost
   :members:

