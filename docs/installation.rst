Installation
============

Funannotate2 can be installed using pip or conda.

Requirements
-----------

Funannotate2 has the following dependencies:

* Python 3.7 or later
* Minimap2
* Miniprot
* SNAP
* Augustus (==3.5.0)
* GlimmerHMM
* GeneMark-ES/ET [optional]
* Diamond
* tRNAscan-SE
* table2asn
* natsort
* numpy
* mappy
* gfftk (>=24.10.29)
* BUSCOlite (>=24.7.29)
* gapmm2
* pyhmmer (>=0.10.15)
* pyfastx (>=2.0.0)
* requests
* gb-io (>=0.3.2)
* json-repair
* pytantan

Some tools like GeneMark-ES/ET must be installed manually due to licensing restrictions (see below).

Installing with conda
------------------

The recommended way to install Funannotate2 is using conda/mamba:

.. code-block:: bash

    mamba create -n funannotate2 gfftk gapmm2 minimap2 miniprot snap "augustus==3.5.0" glimmerhmm diamond trnascan-se table2asn gb-io buscolite
    conda activate funannotate2
    python -m pip install git+https://github.com/nextgenusfs/funannotate2.git


Installation on apple silicon (M series) is a little bit more involved due to some dependency issues and non-native builds of some software.  I've not been able to find or build a version of `augustus` that will run, so instead I've been running `augustus` and `genemark` locally with Docker.  I've setup two repos with instructions on how to get this working (Need Docker Desktop installed) and then will need to put the bash wrapper files in your PATH to mimic the CLI interface.

https://github.com/nextgenusfs/dockerized-augustus

https://github.com/nextgenusfs/dockerized-genemark

Once that is working, you can then install most of the remaining dependencies with conda, although we need to leave out both `buscolite` and `funannotate2` because they have `augustus` as a dependency, instead we will install those python packages with pip.

.. code-block:: bash

    # first install most of the dependencies
    mamba create -n funannotate2 --platform osx-64 "python>=3.7,<3.12" gfftk gapmm2 minimap2 miniprot snap glimmerhmm diamond trnascan-se gb-io pyhmmer pyfastx requests json-repair

    # we can then add the required FUNANNOTATE2_DB env variable to the conda environment, note need to reactivate to use it
    conda activate funannotate2
    conda env config vars set FUNANNOTATE2_DB=/path/to/funannotate2-db
    conda deactivate

    # now reactivate environment, and install the remaining python dependencies with pip
    conda activate funannotate2
    python -m pip install buscolite funannotate2

    # now we can install the databases
    funannotate2 install -d all


Installing with pip
-----------------

To install the latest release version using pip:

.. code-block:: bash

    pip install funannotate2

To install the development version directly from GitHub:

.. code-block:: bash

    pip install git+https://github.com/nextgenusfs/funannotate2.git


Verifying Installation
-------------------

To verify that Funannotate2 is installed correctly:

.. code-block:: bash

    funannotate2 --version

This should display the version of Funannotate2.


Installing Databases
-------------------

Funannotate2 requires several databases to be installed. The funannotate2 scripts expect the $FUNANNOTATE2_DB environment variable to be set. These can be installed using the following command:

.. code-block:: bash

    export FUNANNOTATE2_DB=/path/to/funannotate2_db

    funannotate2 install -d all



Installing GeneMark
----------------

GeneMark-ES/ET must be installed manually due to licensing restrictions:

1. Register and download GeneMark-ES/ET from the `GeneMark website <http://exon.gatech.edu/GeneMark/license_download.cgi>`_
2. Follow the installation instructions provided with the download
3. Make sure the GeneMark executables are in your PATH
