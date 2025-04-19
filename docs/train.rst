Training Ab Initio Gene Predictors
===============================

The ``funannotate2 train`` command trains ab initio gene prediction algorithms to improve gene prediction accuracy. This step is optional but recommended for best results, especially for non-model organisms.

Basic Usage
----------

.. code-block:: bash

    funannotate2 train -f genome.fasta -s "Aspergillus fumigatus" -o train_results

Required Arguments
----------------

* ``-f, --fasta``: Input genome FASTA file
* ``-s, --species``: Species name (e.g., "Aspergillus fumigatus")
* ``-o, --out``: Output folder name

Optional Arguments
----------------

* ``-t, --training-set``: Training set to use in GFF3 format
* ``--strain``: Strain/isolate name
* ``--cpus``: Number of CPUs to use (default: 2)
* ``--optimize-augustus``: Run Augustus mediated optimized training (not recommended) (default: False)
* ``--header-len``: Max length for fasta headers (default: 100)

Training Process
--------------

The ``train`` command performs the following steps:

1. **Prepare the genome**:

   * The genome is loaded and quality checks are performed
   * Headers are checked to ensure they are not too long
   * Contigs are checked for non-IUPAC characters

2. **Generate a training set** (if not provided):

   * BUSCOlite is used to identify conserved single-copy orthologs in the genome
   * The identified orthologs are used to create a training set of gene models
   * The training set is filtered to remove problematic gene models

3. **Split the training set**:

   * The training set is split into test and train sets
   * The test set is used to evaluate the performance of the trained models

4. **Train Augustus**:

   * Augustus is trained using the training set
   * The training process generates species-specific parameters for Augustus
   * The trained parameters are evaluated using the test set

5. **Train SNAP**:

   * SNAP is trained using the training set
   * The training process generates species-specific parameters for SNAP
   * The trained parameters are evaluated using the test set

6. **Train GlimmerHMM**:

   * GlimmerHMM is trained using the training set
   * The training process generates species-specific parameters for GlimmerHMM
   * The trained parameters are evaluated using the test set

7. **Save the trained parameters**:

   * The trained parameters for all tools are saved in a JSON file
   * The JSON file can be used with the ``predict`` command
   * The trained parameters are also saved in the funannotate2 database for future use

Output Files
----------

The ``train`` command generates the following output files in the specified output directory:

* **params.json**: JSON file containing the trained parameters for all tools
* **training-models.train.gff3**: The training set used for training
* **training-models.test.gff3**: The test set used for evaluation

The ``train_misc`` directory contains intermediate files and detailed results from the training process:

* **augustus/**: Directory containing Augustus training files
* **snap/**: Directory containing SNAP training files
* **glimmerhmm/**: Directory containing GlimmerHMM training files
* **busco/**: Directory containing BUSCOlite results (if used to generate the training set)

Using Trained Parameters
---------------------

The trained parameters can be used with the ``predict`` command in two ways:

1. **Using the params.json file**:

   .. code-block:: bash

       funannotate2 predict -f genome.fasta -o predict_results -p train_results/params.json -s "Aspergillus fumigatus"

2. **Using the species name** (if the parameters have been saved in the funannotate2 database):

   .. code-block:: bash

       funannotate2 predict -f genome.fasta -o predict_results -p aspergillus_fumigatus -s "Aspergillus fumigatus"

To save the trained parameters in the funannotate2 database, use the ``species`` command:

.. code-block:: bash

    funannotate2 species -l train_results/params.json

This will make the trained parameters available for future use with the species name as the identifier.
