API Reference
=============

.. toctree::
   :maxdepth: 2

   API/clean
   API/predict
   API/annotate
   API/compare
   API/search
   API/utilities
   API/log

Funannotate2 provides a Python API that can be used to integrate the annotation pipeline into other tools or workflows.

Core Modules
-----------

Funannotate2 consists of several core modules:

- **clean**: Functions for cleaning and preparing genome assemblies
- **predict**: Functions for predicting genes in genome assemblies
- **annotate**: Functions for functionally annotating predicted genes
- **compare**: Functions for comparing multiple genome annotations
- **search**: Functions for searching sequences against various databases
- **utilities**: Utility functions used by other modules
- **log**: Logging functions used by other modules

Using the API
-----------

Here's an example of how to use the Funannotate2 API to clean a genome assembly:

.. code-block:: python

    from funannotate2.clean import clean

    # Create arguments object
    class Args:
        def __init__(self):
            self.input = "genome.fasta"
            self.out = "cleaned_genome.fasta"
            self.minlen = 1000
            self.species = "Aspergillus fumigatus"
            self.strain = "Af293"
            self.header_slice = None
            self.sort = "size"
            self.cpus = 1
            self.tmpdir = None
            self.logfile = None
            self.force = False

    args = Args()

    # Clean the genome
    clean(args)

Here's an example of how to use the Funannotate2 API to predict genes:

.. code-block:: python

    from funannotate2.predict import predict

    # Create arguments object
    class Args:
        def __init__(self):
            self.input = "cleaned_genome.fasta"
            self.out = "predict_results"
            self.species = "Aspergillus fumigatus"
            self.strain = "Af293"
            self.protein_evidence = ["proteins.fasta"]
            self.transcript_evidence = ["transcripts.fasta"]
            self.augustus_species = "aspergillus_fumigatus"
            self.genemark_mode = "ES"
            self.busco_db = "fungi"
            self.busco_seed_species = None
            self.min_intron_len = 10
            self.max_intron_len = 3000
            self.min_protein_len = 50
            self.cpus = 1
            self.tmpdir = None
            self.logfile = None
            self.force = False

    args = Args()

    # Predict genes
    predict(args)

Here's an example of how to use the Funannotate2 API to functionally annotate predicted genes:

.. code-block:: python

    from funannotate2.annotate import annotate

    # Create arguments object
    class Args:
        def __init__(self):
            self.gff3 = "predict_results/funannotate_predict.gff3"
            self.fasta = "cleaned_genome.fasta"
            self.out = "annotate_results"
            self.species = "Aspergillus fumigatus"
            self.strain = "Af293"
            self.pfam = True
            self.dbcan = True
            self.merops = True
            self.swissprot = True
            self.busco = True
            self.busco_db = "fungi"
            self.cpus = 1
            self.tmpdir = None
            self.logfile = None
            self.force = False

    args = Args()

    # Annotate genes
    annotate(args)

Here's an example of how to use the Funannotate2 API to compare multiple genome annotations:

.. code-block:: python

    from funannotate2.compare import compare

    # Create arguments object
    class Args:
        def __init__(self):
            self.input = ["annotate_results", "other_genome_results"]
            self.out = "compare_results"
            self.names = ["Af293", "Other"]
            self.cpus = 1
            self.tmpdir = None
            self.logfile = None
            self.force = False

    args = Args()

    # Compare genomes
    compare(args)
