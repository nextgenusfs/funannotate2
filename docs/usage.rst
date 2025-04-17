Usage
=====

Funannotate2 provides a command-line interface for annotating eukaryotic genomes. The workflow typically consists of the following steps:

1. Clean the genome assembly
2. Predict genes
3. Functionally annotate the predicted genes

Command-Line Interface
---------------------

Funannotate2 provides several commands for different stages of the annotation process:

.. code-block:: bash

    funannotate2 <command> <options>

Available commands:

* ``clean``: Clean and prepare the genome assembly
* ``predict``: Predict genes in the genome
* ``annotate``: Functionally annotate predicted genes
* ``install``: Install and manage databases
* ``train``: Train ab initio gene prediction algorithms
* ``species``: Manage trained species models

Each command has its own set of options. Use ``funannotate2 <command> --help`` to see the available options for each command.

Genome Cleaning
--------------

The first step in the annotation process is to clean and prepare the genome assembly:

.. code-block:: bash

    funannotate2 clean -i genome.fasta -o cleaned_genome.fasta

Options:

* ``-i, --input``: Input genome FASTA file
* ``-o, --out``: Output cleaned genome FASTA file
* ``--minlen``: Minimum contig length to keep (default: 500)
* ``--species``: Species name (e.g., "Aspergillus fumigatus")
* ``--strain``: Strain name (e.g., "Af293")
* ``--header_slice``: Slice contig headers at this character (e.g., " ")
* ``--sort``: Sort contigs by size (default: True)

Gene Prediction
-------------

After cleaning the genome, you can predict genes:

.. code-block:: bash

    funannotate2 predict -i cleaned_genome.fasta -o predict_results -s "Aspergillus fumigatus" --protein_evidence proteins.fasta --transcript_evidence transcripts.fasta

Options:

* ``-i, --input``: Input genome FASTA file
* ``-o, --out``: Output directory
* ``-s, --species``: Species name (e.g., "Aspergillus fumigatus")
* ``--strain``: Strain name (e.g., "Af293")
* ``--protein_evidence``: Protein evidence FASTA file(s)
* ``--transcript_evidence``: Transcript evidence FASTA file(s)
* ``--augustus_species``: Augustus species model to use
* ``--genemark_mode``: GeneMark mode (ES or ET)
* ``--busco_db``: BUSCO database to use for training
* ``--busco_seed_species``: BUSCO seed species for training
* ``--min_intron_len``: Minimum intron length (default: 10)
* ``--max_intron_len``: Maximum intron length (default: 3000)
* ``--min_protein_len``: Minimum protein length (default: 50)
* ``--cpus``: Number of CPUs to use (default: 1)

Functional Annotation
-------------------

After predicting genes, you can functionally annotate them:

.. code-block:: bash

    funannotate2 annotate --gff3 predict_results/funannotate_predict.gff3 --fasta cleaned_genome.fasta -o annotate_results -s "Aspergillus fumigatus"

Options:

* ``--gff3``: Input GFF3 file from predict step
* ``--fasta``: Input genome FASTA file
* ``-o, --out``: Output directory
* ``-s, --species``: Species name (e.g., "Aspergillus fumigatus")
* ``--strain``: Strain name (e.g., "Af293")
* ``--pfam``: Run Pfam annotation (default: True)
* ``--dbcan``: Run dbCAN annotation (default: True)
* ``--merops``: Run MEROPS annotation (default: True)
* ``--swissprot``: Run SwissProt annotation (default: True)
* ``--busco``: Run BUSCO annotation (default: True)
* ``--busco_db``: BUSCO database to use
* ``--cpus``: Number of CPUs to use (default: 1)



Example Workflow
--------------

Here's an example workflow for annotating a fungal genome:

.. code-block:: bash

    # Clean the genome
    funannotate2 clean -i raw_genome.fasta -o cleaned_genome.fasta --minlen 1000 -s "Aspergillus fumigatus" --strain "Af293"

    # Predict genes
    funannotate2 predict -i cleaned_genome.fasta -o predict_results -s "Aspergillus fumigatus" --strain "Af293" \
        --protein_evidence uniprot_fungi.fasta --transcript_evidence rnaseq_transcripts.fasta \
        --augustus_species aspergillus_fumigatus --genemark_mode ES --busco_db fungi --cpus 16

    # Functionally annotate genes
    funannotate2 annotate --gff3 predict_results/funannotate_predict.gff3 --fasta cleaned_genome.fasta \
        -o annotate_results -s "Aspergillus fumigatus" --strain "Af293" \
        --pfam --dbcan --merops --swissprot --busco --busco_db fungi --cpus 16
