Annotate Module
==============

.. automodule:: funannotate2.annotate
   :members:
   :undoc-members:
   :show-inheritance:

The annotate module provides functions for functionally annotating predicted genes.

Key Functions
------------

annotate
~~~~~~~

.. autofunction:: funannotate2.annotate.annotate

This is the main function for functionally annotating predicted genes. It takes a GFF3 file with gene predictions as input and produces various annotation files as output.

check_inputs
~~~~~~~~~~~

.. autofunction:: funannotate2.annotate.check_inputs

Checks that the input files exist and are valid.

find_input_files
~~~~~~~~~~~~~

.. autofunction:: funannotate2.annotate.find_input_files

Finds input files in the specified directory.

parse_annotations
~~~~~~~~~~~~~~

.. autofunction:: funannotate2.annotate.parse_annotations

Parses annotation results from various sources and adds them to the gene models.

naming_slug
~~~~~~~~~

.. autofunction:: funannotate2.annotate.naming_slug

Generates a slug for the species and strain combination.

_sortDict
~~~~~~~

.. autofunction:: funannotate2.annotate._sortDict

Sorts a dictionary by location.
