Search Module
============

.. automodule:: funannotate2.search
   :members:
   :undoc-members:
   :show-inheritance:

The search module provides functions for searching sequences against various databases.

Key Functions
------------

digitize_sequences
~~~~~~~~~~~~~~~

.. autofunction:: funannotate2.search.digitize_sequences

Digitizes sequences for faster searching.

pfam_search
~~~~~~~~~

.. autofunction:: funannotate2.search.pfam_search

Searches protein sequences against the Pfam database.

dbcan_search
~~~~~~~~~~

.. autofunction:: funannotate2.search.dbcan_search

Searches protein sequences against the dbCAN database.

merops_blast
~~~~~~~~~~

.. autofunction:: funannotate2.search.merops_blast

Searches protein sequences against the MEROPS database.

swissprot_blast
~~~~~~~~~~~~

.. autofunction:: funannotate2.search.swissprot_blast

Searches protein sequences against the SwissProt database.

busco_search
~~~~~~~~~~

.. autofunction:: funannotate2.search.busco_search

Searches protein sequences against the BUSCO database.

pfam2tsv
~~~~~~~

.. autofunction:: funannotate2.search.pfam2tsv

Converts Pfam search results to TSV format.

dbcan2tsv
~~~~~~~

.. autofunction:: funannotate2.search.dbcan2tsv

Converts dbCAN search results to TSV format.

merops2tsv
~~~~~~~~

.. autofunction:: funannotate2.search.merops2tsv

Converts MEROPS search results to TSV format.

swissprot2tsv
~~~~~~~~~~~

.. autofunction:: funannotate2.search.swissprot2tsv

Converts SwissProt search results to TSV format.

busco2tsv
~~~~~~~

.. autofunction:: funannotate2.search.busco2tsv

Converts BUSCO search results to TSV format.
