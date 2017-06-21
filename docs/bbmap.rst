.. _bbmap:

******************
Working With BBMap
******************

`BBMap <http://sourceforge.net/projects/bbmap/>`_ is a short read aligner.  We
tested several aligners and found BBMap the most accurate and sensitive.  It's
also relatively fast, and happens to be easy to work with as well.

.. _bbmap_installation:

Installing BBMap
================

First, ensure that you have Java 7 installed.

`Download BBMap <http://sourceforge.net/projects/bbmap/files/latest/download>`_
and extract it to a location of your choice::
    
    # Pick the location for install:
    cd $HOME
    
    # Extract BBMap (extracts to directory 'bbmap')
    tar xzf BBMap_XX.XX.tar.gz

Add it to your PATH::
    
    export PATH=$PATH:$HOME/bbmap

On Mac OS X, you'll also need to set JAVA_HOME::

    export JAVA_HOME=`/usr/libexec/java_home`

Note that BBMap requires Java 7.

.. _building_bbmap_dbs:


Building BBMap Databases
========================

Building a BBMap database is relatively straightforward.  

1. Download a genome or annotation FASTA file.

   * For example, Ensembl provides the `zebrafish genome <ftp://ftp.ensembl.org/pub/release-89/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.toplevel.fa.gz>`_.
   * For example, mirBase provides a `miRNA hairpin annotation FASTA file <ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz>`_.
2. Unzip the zipped FASTA file.  This can sometimes be accomplished by simply double
   clicking on the zipped filename.  This can also be accomplished via
   a terminal window with the following command:

   .. code-block:: bash

      gunzip Danio_rerio.GRCz10.dna.toplevel.fa.gz

3. Create the BBMap database.  Though there are several ways to do this,
   the following approach is known to work with this demo:

   .. code-block:: bash

      bbmap.sh k=7 path=BBMap/Danio_rerio.GRCz10.dna.toplevel \
          ref=Danio_rerio.GRCz10.dna.toplevel.fa.gz

   This may take several minutes depending on your hardware.


BBMap Memory Requirements
=========================

BBMap uses a fair amount of memory.  While it should run fine on any modern
compute node, some laptops and desktops may not have sufficient memory
resources to load an entire reference genome database into BBMap.

