Database Installation
===================

The ``funannotate2 install`` command is used to download and install the databases required for gene prediction and functional annotation.

Basic Usage
----------

To install all databases:

.. code-block:: bash

    funannotate2 install -d all

Required Arguments
----------------

* ``-d, --db``: Databases to install [all,merops,uniprot,dbCAN,pfam,go,mibig,interpro,gene2product,mito]

Optional Arguments
----------------

* ``-s, --show``: Show currently installed databases (default: False)
* ``-w, --wget``: Use wget for downloading (default: False)
* ``-f, --force``: Force re-download/re-install of all databases (default: False)
* ``-u, --update``: Update databases if change detected (default: False)

Available Databases
----------------

The following databases can be installed:

* **all**: Install all databases
* **merops**: MEROPS protease database
* **uniprot**: UniProtKB/Swiss-Prot database
* **dbCAN**: Database for Carbohydrate-active enzyme ANnotation
* **pfam**: Pfam protein families database
* **go**: Gene Ontology database
* **mibig**: Minimum Information about a Biosynthetic Gene cluster
* **interpro**: InterPro protein families database
* **gene2product**: Curated gene names and product descriptions
* **mito**: Mitochondrial genome database

Showing Installed Databases
------------------------

To see which databases are currently installed:

.. code-block:: bash

    funannotate2 install -s

This will display a list of installed databases, their versions, and installation dates.

Updating Databases
---------------

To update all installed databases:

.. code-block:: bash

    funannotate2 install -d all -u

This will check each database for updates and download new versions if available.

Forcing Reinstallation
-------------------

To force reinstallation of a database:

.. code-block:: bash

    funannotate2 install -d pfam -f

This will delete the existing database and download it again, regardless of whether it's already installed.

Database Location
--------------

Databases are installed in the directory specified by the ``$FUNANNOTATE2_DB`` environment variable. This variable must be set before running the ``install`` command.

To set the environment variable:

.. code-block:: bash

    # For bash/zsh
    export FUNANNOTATE2_DB=/path/to/database/directory

    # For csh/tcsh
    setenv FUNANNOTATE2_DB /path/to/database/directory

You can add this line to your shell's startup file (e.g., ~/.bashrc, ~/.zshrc) to make it permanent.

Database Information
-----------------

Information about installed databases is stored in the ``funannotate-db-info.json`` file in the database directory. This file contains:

* Database names
* Installation dates
* Version information
* File paths

This information is used by other funannotate2 commands to locate the required databases.

Troubleshooting
------------

If you encounter issues with database installation:

1. **Check disk space**: Ensure you have enough disk space for the databases (several GB may be required)
2. **Check permissions**: Ensure you have write permissions to the database directory
3. **Check internet connection**: Ensure you have a stable internet connection
4. **Try wget**: Use the ``-w`` option to use wget instead of the default downloader
5. **Check logs**: Look for error messages in the output

If problems persist, try installing databases individually rather than all at once.
