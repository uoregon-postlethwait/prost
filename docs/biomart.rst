.. _biomart:

*********************************************************
Generating a BioMart Other RNA annotations BBMap database
*********************************************************

This page describes how to generate an annotations database containing "Other
RNAs" such as snoRNAs and lncRNAs.

If you do not wish to supply a real BioMart file, simply download and point 
*Prost!* at the following fake BioMart file: `fake_biomart_file.fa <https://raw.githubusercontent.com/uoregon-postlethwait/prost/master/fake_biomart_file.fa>`_

Download the sequences from BioMart
```````````````````````````````````

+ Open http://ensembl.org/biomart in your browser.
+ Make the following selections:

  + Click **- CHOOSE DATABASE -** and select **Ensembl Genes xx** where 
    *xx* is the current version of Ensembl.
  + Click **- CHOOSE DATASET -** and select your species of interest (e.g.
    *Danio rerio genes (GRCz10)*).
  + Click **Filters**

    - Click the **+** button next to **GENE** to expand that box.
    - Under **Gene Type**, select all of the following if present (you may have
      to ctrl-click or cmd-click to select multiples at the same time):

      * lincRNA
      * miRNA
      * misc_RNA
      * Mt_rRNA
      * Mt_tRNA
      * rRNA
      * snoRNA
      * snRNA

  - Click **Attributes** 

    - Click the **Sequences** radio button.
    - Click the **+** button next to **SEQUENCES** to expand that box.

      - Click the **cDNA sequences** radio button.

    - Click the **+** button next to **HEADER INFORMATION** to expand that box.

      - Uncheck all the checkboxes.
      - Check the these boxes **in the following order**: 

        - Gene Name
        - Transcript type
        - Transcript stable ID

+ *Optional*: Click *Count* to see how many sequences will be returned.
+ Click *Results* to generate the fasta file.
+ Click *Go* to download the fasta file containing the other RNAs.


- Move the downloaded file (typically named ``mart_export.txt``) somewhere
  appropriate on your computer (or onto another computer if that is where you'll
  be running *Prost!*). 
- Rename the downloaded file (typically named ``mart_export.txt``) to something
  more appropriate (e.g. ``YOUR_SPECIES_biomart_other_rnas.fa``).
- Build a BBMap database from the FASTA file (see :ref:`building_bbmap_dbs`).
- Update the *Prost!* configuration file (usually ``prost.config``) to reflect
  the newly created BBMap database, for example:

.. code-block:: ini

   [AnnotationAlignment3]
   type: BiomartOtherRNAAnnotation
   name: other
   tool: bbmap
   db: /path/to/YOUR_SPECIES_biomart_other_rnas
   max_3p_mismatches: 0
   max_non_3p_mismatches: 0
   allow_indels: no
