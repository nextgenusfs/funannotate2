Managing Trained Species
=====================

The ``funannotate2 species`` command allows you to view and manage trained species models in the funannotate2 database. These trained models are used by the ``predict`` command to improve gene prediction accuracy.

Basic Usage
----------

To list all trained species in the database:

.. code-block:: bash

    funannotate2 species

Optional Arguments
----------------

* ``-l, --load``: Load a new species with a *.params.json file
* ``-d, --delete``: Delete a species from database
* ``-f, --format``: Format to show existing species in (default: table)

Viewing Trained Species
--------------------

By default, the ``species`` command lists all trained species in the database in a tabular format:

.. code-block:: bash

    funannotate2 species

You can change the output format to JSON or YAML:

.. code-block:: bash

    funannotate2 species -f json
    funannotate2 species -f yaml

Loading a New Species
------------------

After training a species using the ``train`` command, you can add it to the database:

.. code-block:: bash

    funannotate2 species -l /path/to/params.json

The params.json file is typically found in the output directory of the ``train`` command:

.. code-block:: bash

    funannotate2 species -l train_results/params.json

Deleting a Species
---------------

You can remove a species from the database:

.. code-block:: bash

    funannotate2 species -d species_name

For example:

.. code-block:: bash

    funannotate2 species -d aspergillus_fumigatus

How Species Models are Used
------------------------

When you run the ``predict`` command, you can specify a trained species model:

.. code-block:: bash

    funannotate2 predict -f genome.fasta -o predict_results -p species_name -s "Species name"

For example:

.. code-block:: bash

    funannotate2 predict -f genome.fasta -o predict_results -p aspergillus_fumigatus -s "Aspergillus fumigatus"

The trained species model provides parameters for ab initio gene predictors:

1. **Augustus**: Species-specific parameters for splice sites, start/stop codons, etc.
2. **SNAP**: Species-specific HMM parameters
3. **GlimmerHMM**: Species-specific parameters for gene structure

Using a species model that is closely related to your target organism can significantly improve gene prediction accuracy.

Pretrained Species
---------------

Funannotate2 comes with several pretrained species models for common organisms. You can see the list of available pretrained species with the ``species`` command.

If your organism is not in the list, you can:

1. Use a closely related species model
2. Train a new model using the ``train`` command
3. Load the new model into the database using the ``species -l`` command

Species Model Storage
------------------

Species models are stored in the funannotate2 database directory, which is specified by the ``$FUNANNOTATE2_DB`` environment variable. Each species model includes:

* Parameters for Augustus
* Parameters for SNAP
* Parameters for GlimmerHMM
* Metadata about the training process

The models are stored in a structured format that allows them to be easily loaded and used by the ``predict`` command.
