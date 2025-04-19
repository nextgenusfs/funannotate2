Genome Cleaning
==============

The ``funannotate2 clean`` command prepares a genome assembly for annotation by removing duplicated contigs, sorting contigs by size, and optionally renaming contig headers.

Basic Usage
----------

.. code-block:: bash

    funannotate2 clean -f genome.fasta -o cleaned_genome.fasta

Required Arguments
----------------

* ``-f, --fasta``: Input genome FASTA file
* ``-o, --out``: Output cleaned genome FASTA file

Optional Arguments
----------------

* ``-p, --pident``: Percent identity threshold for identifying duplicated contigs (default: 95)
* ``-c, --cov``: Coverage threshold for identifying duplicated contigs (default: 95)
* ``-m, --minlen``: Minimum contig length to keep (default: 500)
* ``-r, --rename``: Rename contigs with this basename (e.g., "scaffold_")
* ``--cpus``: Number of CPUs to use (default: 2)
* ``--tmpdir``: Directory for temporary files
* ``--exhaustive``: Compute every contig, else stop at N50 (default: False)
* ``--logfile``: Write logs to file
* ``--debug``: Debug the output (default: False)

Cleaning Process
--------------

The ``clean`` command performs the following steps:

1. **Load and sort the genome**:

   * The genome is loaded from the input FASTA file
   * Contigs are sorted by length, from shortest to longest
   * Basic statistics are calculated, including N50

2. **Filter contigs by length**:

   * Contigs shorter than the minimum length are removed
   * This helps eliminate small, potentially problematic contigs

3. **Check for duplicated contigs**:

   * Starting with the smallest contigs, each contig is aligned against all larger contigs using `minimap2 <https://github.com/lh3/minimap2>`_
   * If a contig is found to be duplicated elsewhere in the genome (based on percent identity and coverage thresholds), it is marked for removal
   * By default, the process stops at the N50 contig size to save time, as larger contigs are less likely to be duplicated
   * If the ``--exhaustive`` option is used, all contigs are checked for duplication

4. **Rename contigs** (if requested):

   * If the ``--rename`` option is used, contigs are renamed with the provided basename followed by a number (e.g., "scaffold_1", "scaffold_2", etc.)
   * Contigs are numbered from largest to smallest

5. **Write the cleaned genome**:

   * The cleaned genome is written to the output file
   * Duplicated contigs are excluded
   * Contigs are written in order from largest to smallest

Why Clean Your Genome?
--------------------

Cleaning your genome assembly before annotation offers several benefits:

1. **Removes redundant sequences**:

   * Duplicated contigs can lead to redundant gene annotations
   * Removing duplicates ensures each gene is annotated only once

2. **Improves annotation quality**:

   * Small, fragmented contigs often contain partial genes or repetitive elements
   * Removing these contigs can improve the overall quality of gene predictions

3. **Standardizes contig names**:

   * Consistent, simple contig names make downstream analysis easier
   * Some tools have limitations on header length or format

4. **Reduces computational requirements**:

   * Fewer contigs means faster processing in subsequent steps
   * Removing small contigs can significantly reduce the number of contigs without losing much sequence

Example Usage
-----------

Basic cleaning with default parameters:

.. code-block:: bash

    funannotate2 clean -f raw_genome.fasta -o cleaned_genome.fasta

Cleaning with custom parameters:

.. code-block:: bash

    funannotate2 clean -f raw_genome.fasta -o cleaned_genome.fasta -m 1000 -p 98 -c 98 -r scaffold_

Exhaustive cleaning (check all contigs for duplication):

.. code-block:: bash

    funannotate2 clean -f raw_genome.fasta -o cleaned_genome.fasta --exhaustive

Output
-----

The ``clean`` command produces a single output file: the cleaned genome in FASTA format. The command also outputs statistics about the cleaning process, including:

* Number of contigs in the input genome
* Number of contigs larger than the minimum length
* N50 of the input genome
* Number of duplicated contigs found
* Number of contigs written to the output file
