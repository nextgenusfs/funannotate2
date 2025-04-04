Installation
============

Funannotate2 can be installed using pip or conda.

Requirements
-----------

Funannotate2 has the following dependencies:

* Python 3.8 or later
* Biopython
* EvidenceModeler (>=2)
* Minimap2
* Miniprot
* SNAP
* Augustus (==3.5.0)
* GlimmerHMM
* Diamond
* tRNAscan-SE
* table2asn

Some tools like GeneMark-ES/ET must be installed manually due to licensing restrictions.

Installing with conda
------------------

The recommended way to install Funannotate2 is using conda/mamba:

.. code-block:: bash

    mamba create -n funannotate2 "python<=3.10" biopython "evidencemodeler>=2" minimap2 miniprot snap "augustus==3.5.0" glimmerhmm diamond trnascan-se table2asn
    conda activate funannotate2
    python -m pip install git+https://github.com/nextgenusfs/funannotate2.git

Installing with pip
-----------------

To install the latest release version using pip:

.. code-block:: bash

    pip install funannotate2

To install the development version directly from GitHub:

.. code-block:: bash

    pip install git+https://github.com/nextgenusfs/funannotate2.git

Installing GeneMark
----------------

GeneMark-ES/ET must be installed manually due to licensing restrictions:

1. Register and download GeneMark-ES/ET from the `GeneMark website <http://exon.gatech.edu/GeneMark/license_download.cgi>`_
2. Follow the installation instructions provided with the download
3. Make sure the GeneMark executables are in your PATH

Verifying Installation
-------------------

To verify that Funannotate2 is installed correctly:

.. code-block:: bash

    funannotate2 --version

This should display the version of Funannotate2.

To check if all dependencies are installed correctly:

.. code-block:: bash

    funannotate2 check --dependencies

This will check if all required dependencies are installed and available in your PATH.
