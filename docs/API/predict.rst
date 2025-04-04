Predict Module
=============

.. automodule:: funannotate2.predict
   :members:
   :undoc-members:
   :show-inheritance:

The predict module provides functions for predicting genes in genome assemblies.

Key Functions
------------

predict
~~~~~~

.. autofunction:: funannotate2.predict.predict

This is the main function for predicting genes in genome assemblies. It takes a genome FASTA file as input and produces a GFF3 file with gene predictions as output.

check_inputs
~~~~~~~~~~~

.. autofunction:: funannotate2.predict.check_inputs

Checks that the input files exist and are valid.

run_genemark_es
~~~~~~~~~~~~~

.. autofunction:: funannotate2.predict.run_genemark_es

Runs GeneMark-ES on the genome assembly.

run_augustus
~~~~~~~~~~

.. autofunction:: funannotate2.predict.run_augustus

Runs Augustus on the genome assembly.

run_miniprot
~~~~~~~~~~

.. autofunction:: funannotate2.predict.run_miniprot

Runs Miniprot to align protein evidence to the genome assembly.

run_minimap2_transcripts
~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: funannotate2.predict.run_minimap2_transcripts

Runs Minimap2 to align transcript evidence to the genome assembly.

run_minimap2_proteins
~~~~~~~~~~~~~~~~~~

.. autofunction:: funannotate2.predict.run_minimap2_proteins

Runs Minimap2 to align protein evidence to the genome assembly.

merge_predictions
~~~~~~~~~~~~~~

.. autofunction:: funannotate2.predict.merge_predictions

Merges gene predictions from multiple sources.
