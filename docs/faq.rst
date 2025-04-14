Frequently Asked Questions
==========================

General Questions
---------------

What is Funannotate2?
~~~~~~~~~~~~~~~~~~~

Funannotate2 is a comprehensive eukaryotic genome annotation pipeline that provides a complete workflow for annotating fungal, plant, and other eukaryotic genomes. It integrates various tools and databases to produce high-quality gene predictions and functional annotations.

What types of genomes can Funannotate2 annotate?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Funannotate2 is designed primarily for fungal genomes, but it can also be used to annotate other eukaryotic genomes such as plants, insects, and other organisms. However, the default parameters and databases are optimized for fungal genomes.

What are the system requirements for Funannotate2?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Funannotate2 requires:

- Linux or macOS operating system
- Python 3.7 or later
- At least 8 GB of RAM (16 GB or more recommended for larger genomes)
- At least 50 GB of free disk space
- Multiple CPU cores (8 or more recommended for faster processing)

Installation Questions
-------------------

How do I install Funannotate2?
~~~~~~~~~~~~~~~~~~~~~~~~~~~

See the :doc:`installation` page for detailed installation instructions.

Why can't I install GeneMark-ES through conda?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GeneMark-ES is not available through conda due to licensing restrictions. You need to register and download it manually from the `GeneMark website <http://exon.gatech.edu/GeneMark/license_download.cgi>`_.


Usage Questions
------------

What is the recommended workflow for annotating a genome?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The recommended workflow is:

1. Clean the genome assembly using ``funannotate2 clean``
2. Train ab initio prediction tools using ``funannotate2 train``
3. Predict genes using ``funannotate2 predict``
4. Functionally annotate the predicted genes using ``funannotate2 annotate``

See the :doc:`tutorial` for a detailed example.

How can I improve gene prediction accuracy?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To improve gene prediction accuracy:

1. Use high-quality protein evidence from closely related species
2. Use transcript evidence from RNA-seq data
3. Use a species-specific Augustus model
4. Use the appropriate GeneMark mode (ES for self-training, ET for transcript-guided)
5. Use a BUSCO database appropriate for your organism

What databases does Funannotate2 use for functional annotation?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Funannotate2 uses the following databases for functional annotation:

- Pfam: Protein domain annotations
- dbCAN: Carbohydrate-active enzyme annotations
- MEROPS: Peptidase annotations
- SwissProt: Protein annotations
- BUSCO: Benchmarking Universal Single-Copy Orthologs

How can I add custom functional annotations?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can add custom functional annotations by:

1. Creating a custom database in the appropriate format (FASTA, HMM, etc.)
2. Using the appropriate search tool (BLAST, HMMER, etc.) to search your proteins against the custom database
3. Parsing the search results and adding the annotations to the gene models
4. Using the Funannotate2 API to integrate the custom annotations into the annotation pipeline

Troubleshooting
-------------

Why does GeneMark-ES fail on my genome?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GeneMark-ES may fail for several reasons:

1. The genome assembly is too fragmented (try filtering out short contigs)
2. The genome assembly contains too many Ns (try cleaning the genome)
3. The genome is not from a eukaryotic organism (GeneMark-ES is designed for eukaryotes)
4. GeneMark-ES is not installed correctly (check the installation)

Why does Augustus fail on my genome?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Augustus may fail for several reasons:

1. The species model does not exist (try using a different species model)
2. The species model is not appropriate for your organism (try using a more closely related species)
3. Augustus is not installed correctly (check the installation)
4. The genome assembly is too fragmented (try filtering out short contigs)

Why are some of my gene models incomplete?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Gene models may be incomplete for several reasons:

1. The genome assembly is fragmented, and genes span contig boundaries
2. The gene prediction tools failed to identify the complete gene structure
3. The gene is genuinely partial (e.g., pseudogene)

Try using protein and transcript evidence to improve gene model completeness.

How can I report a bug or request a feature?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can report bugs or request features by opening an issue on the `GitHub repository <https://github.com/nextgenusfs/funannotate2/issues>`_.
