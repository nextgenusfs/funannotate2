Usage
=====

Funannotate2 provides a command-line interface for annotating eukaryotic genomes. The workflow typically consists of the following steps:

1. Clean the genome assembly
2. Train ab initio gene prediction algorithms (optional)
3. Predict genes
4. Functionally annotate the predicted genes

Command-Line Interface
---------------------

Funannotate2 provides several commands for different stages of the annotation process:

.. code-block:: bash

    funannotate2 <command> <options>

Available commands:

* ``clean``: Find and remove duplicated contigs, sort by size, rename headers
* ``train``: Train ab initio gene prediction algorithms
* ``predict``: Predict primary gene models in a eukaryotic genome
* ``annotate``: Add functional annotation to gene models
* ``install``: Install and manage databases
* ``species``: View and manage trained species models

Each command has its own set of options. Use ``funannotate2 <command> --help`` to see the available options for each command.

Genome Cleaning
--------------

The first step in the annotation process is to clean and prepare the genome assembly. The ``clean`` command sorts contigs by size, identifies and removes duplicated contigs, and optionally renames contig headers.

.. code-block:: bash

    funannotate2 clean -f genome.fasta -o cleaned_genome.fasta

Required arguments:

* ``-f, --fasta``: Input genome FASTA file
* ``-o, --out``: Output cleaned genome FASTA file

Optional arguments:

* ``-p, --pident``: Percent identity threshold for identifying duplicated contigs (default: 95)
* ``-c, --cov``: Coverage threshold for identifying duplicated contigs (default: 95)
* ``-m, --minlen``: Minimum contig length to keep (default: 500)
* ``-r, --rename``: Rename contigs with this basename (e.g., "scaffold_")
* ``--cpus``: Number of CPUs to use (default: 2)
* ``--tmpdir``: Directory for temporary files
* ``--exhaustive``: Compute every contig, else stop at N50 (default: False)
* ``--logfile``: Write logs to file
* ``--debug``: Debug the output (default: False)

**How it works:**

1. The genome is loaded and contigs are sorted by length
2. Contigs smaller than the minimum length are filtered out
3. Starting with the smallest contigs, each contig is aligned against all larger contigs using minimap2
4. If a contig is found to be duplicated elsewhere in the genome (based on percent identity and coverage thresholds), it is marked for removal
5. By default, the process stops at the N50 contig size to save time, as larger contigs are less likely to be duplicated
6. If the ``--exhaustive`` option is used, all contigs are checked for duplication
7. If the ``--rename`` option is used, contigs are renamed with the provided basename followed by a number (e.g., "scaffold_1", "scaffold_2", etc.)
8. The cleaned genome is written to the output file

Training Ab Initio Gene Predictors
-----------------------------

Before predicting genes, you can train ab initio gene prediction algorithms to improve accuracy. This step is optional but recommended for best results, especially for non-model organisms. You can re-use training data by passing a pretrained species slug or a params.json file.

.. code-block:: bash

    funannotate2 train -f cleaned_genome.fasta -s "Aspergillus nidulans" -o anid_f2

Required arguments:

* ``-f, --fasta``: Input genome FASTA file
* ``-s, --species``: Species name (e.g., "Aspergillus fumigatus")
* ``-o, --out``: Output folder name

Optional arguments:

* ``-t, --training-set``: Training set to use in GFF3 format
* ``--strain``: Strain/isolate name
* ``--cpus``: Number of CPUs to use (default: 2)
* ``--optimize-augustus``: Run Augustus mediated optimized training (not recommended) (default: False)
* ``--header-len``: Max length for fasta headers (default: 100)

**How it works:**

1. If no training set is provided, BUSCOlite is used to identify conserved single-copy orthologs in the genome
2. The identified orthologs are used to create a training set of gene models
3. The training set is filtered and split into test and train sets
4. Augustus, SNAP, and GlimmerHMM are trained using the training set
5. The trained parameters are saved in a JSON file that can be used with the ``predict`` command
6. The trained parameters are also saved in the funannotate2 database for future use

Gene Prediction
-------------

After cleaning the genome and optionally training ab initio gene predictors, you can predict genes:

.. code-block:: bash

    funannotate2 predict -i anid_f2

Required arguments:

* ``-i, --input-dir``: funannotate2 output directory
* ``-f, --fasta``: Input genome FASTA file (softmasked repeats)
* ``-o, --out``: Output folder name
* ``-p, --params, --pretrained``: Params.json or pretrained species slug. Use ``funannotate2 species`` to see pretrained species
* ``-s, --species``: Species name (e.g., "Aspergillus fumigatus")

Optional arguments:

* ``-st, --strain``: Strain/isolate name (e.g., "Af293")
* ``-e, --external``: External gene models/annotation in GFF3 format
* ``-w, --weights``: Gene predictors and weights
* ``-ps, --proteins``: Proteins to use for evidence
* ``-ts, --transcripts``: Transcripts to use for evidence
* ``-c, --cpus``: Number of CPUs to use (default: 2)
* ``-mi, --max-intron``: Maximum intron length (default: 3000)
* ``-hl, --header-len``: Max length for fasta headers (default: 100)
* ``-l, --locus-tag``: Locus tag for genes, perhaps assigned by NCBI, e.g. VC83 (default: FUN2_)
* ``-n, --numbering``: Specify start of gene numbering (default: 1)
* ``--tmpdir``: Volume to write tmp files (default: /tmp)

**How it works:**

1. The genome is analyzed for assembly statistics and softmasked regions
2. If the genome is not softmasked, pytantan is used to quickly softmask repeats
3. If protein evidence is provided, it is aligned to the genome using ``miniprot``
4. If transcript evidence is provided, it is aligned to the genome using ``gapmm2``
5. Evidence alignments are converted to hints for ``augustus``
6. Ab initio gene predictors (Augustus, GeneMark [optional], SNAP, GlimmerHMM) are run on the genome
7. tRNAscan-SE is run to identify tRNA genes
8. The GFFtk consensus module is used to generate consensus gene models from all evidence and ab initio predictions
9. The consensus gene models are filtered and annotated
10. The final gene models are output in GFF3, TBL, and GenBank formats
11. Protein and transcript sequences are extracted from the gene models
12. Summary statistics are generated for the annotation

Functional Annotation
-------------------

After predicting genes, you can functionally annotate them:

.. code-block:: bash

    funannotate2 annotate -i anid_f2

Required arguments:

* ``-i, --input-dir``: funannotate2 output directory
* ``-f, --fasta``: Genome in FASTA format (required if not using --input-dir)
* ``-t, --tbl``: Genome annotation in TBL format (required if not using --input-dir and not using --gff3)
* ``-g, --gff3``: Genome annotation in GFF3 format (required if not using --input-dir and not using --tbl)
* ``-o, --out``: Output folder name (required if not using --input-dir)

Optional arguments:

* ``-a, --annotations``: Annotations files, 3 column TSV [transcript-id, feature, data]
* ``-s, --species``: Species name (e.g., "Aspergillus fumigatus")
* ``-st, --strain``: Strain/isolate name
* ``--cpus``: Number of CPUs to use (default: 2)
* ``--tmpdir``: Volume to write tmp files (default: /tmp)
* ``--curated-names``: Path to custom file with gene-specific annotations (tab-delimited: gene_id annotation_type annotation_value)

**How it works:**

1. The gene models are loaded from the input directory or specified files
2. Protein sequences are extracted from the gene models
3. The proteins are searched against various databases for functional annotation:
   - Pfam-A database using pyhmmer for protein domains
   - dbCAN database using pyhmmer for carbohydrate-active enzymes (CAZymes)
   - UniProtKB/Swiss-Prot database using diamond for protein function
   - MEROPS database using diamond for proteases
   - BUSCOlite for conserved orthologs
4. Gene names and product descriptions are cleaned using a curated database
5. Custom annotations can be provided to override automatic cleaning
6. The annotations are merged into the gene models
7. The annotated gene models are output in GFF3, TBL, and GenBank formats
8. Protein and transcript sequences are extracted from the annotated gene models
9. Summary statistics are generated for the annotation

For more details on using custom curated gene names and products, see the :doc:`annotate` page.



Example Workflow
--------------

Here's an example workflow for annotating a fungal genome:

.. code-block:: bash

    # Clean the genome
    funannotate2 clean -f raw_genome.fasta -o cleaned_genome.fasta -m 1000 -r scaffold_

    # Train ab initio gene predictors (optional)
    funannotate2 train -f cleaned_genome.fasta -s "Aspergillus fumigatus" -o f2_output --strain "Af293" --cpus 16

    # Predict genes using trained parameters
    funannotate2 predict -i f2_output -ps uniprot_fungi.fasta -ts rnaseq_transcripts.fasta --cpus 16

    # Or predict genes using pretrained species
    funannotate2 predict -i f2_output -p aspergillus_fumigatus -ps uniprot_fungi.fasta -ts rnaseq_transcripts.fasta --cpus 16

    # Functionally annotate genes
    funannotate2 annotate -i f2_output --cpus 16

    # Add custom gene/product annotations (optional)
    funannotate2 annotate -i f2_output --cpus 16 --curated-names custom_annotations.txt

For more detailed examples and explanations, see the :doc:`tutorial` page.

Database Installation
------------------

Funannotate2 requires several databases for gene prediction and functional annotation. You can install these databases using the ``install`` command:

.. code-block:: bash

    funannotate2 install -d all

Required arguments:

* ``-d, --db``: Databases to install [all,merops,uniprot,dbCAN,pfam,go,mibig,interpro,gene2product,mito]

Optional arguments:

* ``-s, --show``: Show currently installed databases (default: False)
* ``-w, --wget``: Use wget for downloading (default: False)
* ``-f, --force``: Force re-download/re-install of all databases (default: False)
* ``-u, --update``: Update databases if change detected (default: False)

**How it works:**

1. The command checks for the ``$FUNANNOTATE2_DB`` environment variable, which should point to the directory where databases will be installed
2. If the specified databases are already installed, the command will skip them unless ``--force`` or ``--update`` is used
3. The databases are downloaded from their respective sources and processed for use with funannotate2
4. A record of installed databases is kept in the ``funannotate-db-info.json`` file in the database directory

Managing Trained Species
---------------------

Funannotate2 maintains a database of trained species parameters for gene prediction. You can view and manage these species using the ``species`` command:

.. code-block:: bash

    funannotate2 species

Optional arguments:

* ``-l, --load``: Load a new species with a *.params.json file
* ``-d, --delete``: Delete a species from database
* ``-f, --format``: Format to show existing species in (default: table)

**How it works:**

1. Without arguments, the command lists all trained species in the database
2. With ``--load``, the command adds a new species to the database from a params.json file (typically generated by the ``train`` command)
3. With ``--delete``, the command removes a species from the database
4. The ``--format`` option controls how the species list is displayed (table, json, or yaml)
