Clean Module
===========

.. automodule:: funannotate2.clean
   :members:
   :undoc-members:
   :show-inheritance:

The clean module provides functions for cleaning and preparing genome assemblies for annotation.

Key Functions
------------

clean
~~~~~

.. autofunction:: funannotate2.clean.clean

This is the main function for cleaning and preparing genome assemblies. It takes a genome FASTA file as input and produces a cleaned FASTA file as output.

check_inputs
~~~~~~~~~~~

.. autofunction:: funannotate2.clean.check_inputs

Checks that the input files exist and are valid.

clean_header
~~~~~~~~~~

.. autofunction:: funannotate2.clean.clean_header

Cleans FASTA headers by removing unwanted characters and optionally slicing at a specified character.

filter_contigs
~~~~~~~~~~~~

.. autofunction:: funannotate2.clean.filter_contigs

Filters contigs by minimum length.

sort_contigs
~~~~~~~~~~

.. autofunction:: funannotate2.clean.sort_contigs

Sorts contigs by size or name.
