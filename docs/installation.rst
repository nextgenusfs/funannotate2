Installation
============

Funannotate2 can be installed with Docker (recommended), with `pixi <https://pixi.sh/>`_, or with conda/pip.

Using the Docker image
----------------------

A pre-built image containing ``funannotate2``, ``funannotate2-addons``, ``helixerlite``, all bioconda tooling, and the databases from ``funannotate2 install -d all`` is published on each tagged release to:

* Docker Hub: ``nextgenusfs/funannotate2``
* GHCR: ``ghcr.io/nextgenusfs/funannotate2``

Tags follow SemVer (e.g. ``:26.2.12``, ``:26.2``) with ``:latest`` always pointing at the most recent release.

.. code-block:: bash

    # pull the latest image (~8 GB; databases are baked in)
    docker pull nextgenusfs/funannotate2:latest

    # sanity checks
    docker run --rm nextgenusfs/funannotate2:latest funannotate2 --version
    docker run --rm nextgenusfs/funannotate2:latest funannotate2 install -s

    # run against a local data directory, persisting the BUSCO cache across runs
    mkdir -p $PWD/data $PWD/busco_cache
    docker run --rm -it \
        -v $PWD/data:/data \
        -v $PWD/busco_cache:/opt/busco_cache \
        -e BUSCO_DOWNLOAD_PATH=/opt/busco_cache \
        nextgenusfs/funannotate2:latest \
        funannotate2 predict -i /data/genome.fa -o /data/out --species "My species"

Notes:

* The image is ``linux/amd64`` only. On Apple Silicon it runs under Rosetta 2 emulation (Docker Desktop handles this automatically).
* BUSCO lineages are **not** bundled — at ~90 GB uncompressed they exceed Docker Hub's size ceiling. Mount a host directory (as shown above) so lineages download once and are reused.
* ``GeneMark`` is not included due to licensing. If you need it, install it on the host and mount it into the container.
* ``FUNANNOTATE2_DB`` is already set to ``/opt/funannotate2_db`` inside the image; no host setup required unless you want to override it.

Installing with pixi
--------------------

If you prefer a native install on Linux without Docker, the repository ships a `pixi <https://pixi.sh/>`_ workspace (``pixi.toml`` / ``pixi.lock``) that resolves the same environment used to build the Docker image.

.. code-block:: bash

    # install pixi once (see https://pixi.sh/latest/#installation for alternatives)
    curl -fsSL https://pixi.sh/install.sh | bash

    # clone the repo and install the locked environment
    git clone https://github.com/nextgenusfs/funannotate2.git
    cd funannotate2
    pixi install --locked

    # activate the environment
    pixi shell

    # set the database location and install databases
    export FUNANNOTATE2_DB=/path/to/funannotate2-db
    funannotate2 install -d all

The pixi environment is currently defined for ``linux-64`` only. macOS users should use the Docker image.

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
* pyhmmer (>=0.12.0)
* pyfastx (>=2.0.0)
* requests
* gb-io (>=0.3.2)
* json-repair
* pytantan

Some tools like GeneMark-ES/ET must be installed manually due to licensing restrictions (see below).

Installing with conda
------------------

Funannotate2 can also be installed with conda/mamba:

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
    mamba create -n funannotate2 --platform osx-64 "python>=3.7,<3.13" gfftk gapmm2 minimap2 miniprot snap glimmerhmm diamond trnascan-se gb-io pyhmmer pyfastx requests json-repair "mkl<2022" pytantan

    # we can then add the required FUNANNOTATE2_DB env variable to the conda environment, note need to reactivate to use it
    conda activate funannotate2
    conda env config vars set FUNANNOTATE2_DB=/path/to/funannotate2-db
    conda env config vars set AUGUSTUS_CONFIG_PATH=/path/to/augustus-3.5.0/config
    conda deactivate

    # now reactivate environment, and install the remaining python dependencies with pip
    conda activate funannotate2
    python -m pip install buscolite git+https://github.com/nextgenusfs/funannotate2.git

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

Funannotate2 requires several databases to be installed. Note: funannotate2 scripts expect the $FUNANNOTATE2_DB environment variable to be set. These can be installed using the following command:

.. code-block:: bash

    funannotate2 install -d all



Installing GeneMark
----------------

GeneMark-ES/ET must be installed manually due to licensing restrictions:

1. Register and download GeneMark-ES/ET from the `GeneMark website <http://exon.gatech.edu/GeneMark/license_download.cgi>`_
2. Follow the installation instructions provided with the download
3. Make sure the GeneMark executables are in your PATH
4. You may also need to install GeneMark specific perl libraries, specifically ``perl-hash-merge`` and ``perl-mce`` have been mentioned by users adding to the existing conda environment.
