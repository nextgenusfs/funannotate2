Tutorial
========

This tutorial will guide you through the process of annotating a fungal genome using Funannotate2.

Prerequisites
------------

Before starting this tutorial, make sure you have:

1. Installed Funannotate2 and its dependencies (see :doc:`installation`)
2. A genome assembly in FASTA format
3. Protein evidence (optional but recommended)
4. Transcript evidence (optional but recommended)

Step 1: Clean the Genome Assembly
--------------------------------

The first step is to clean and prepare the genome assembly for annotation:

.. code-block:: bash

    funannotate2 clean -i raw_genome.fasta -o cleaned_genome.fasta --minlen 1000 -s "Aspergillus fumigatus" --strain "Af293"

This command will:

1. Remove contigs shorter than 1000 bp
2. Clean contig headers
3. Sort contigs by size
4. Output the cleaned genome to ``cleaned_genome.fasta``

Step 2: Predict Genes
-------------------

The next step is to predict genes in the cleaned genome assembly:

.. code-block:: bash

    funannotate2 predict -i cleaned_genome.fasta -o predict_results -s "Aspergillus fumigatus" --strain "Af293" \
        --protein_evidence uniprot_fungi.fasta --transcript_evidence rnaseq_transcripts.fasta \
        --augustus_species aspergillus_fumigatus --genemark_mode ES --busco_db fungi --cpus 16

This command will:

1. Run GeneMark-ES to predict genes
2. Run Augustus to predict genes
3. Align protein evidence using Miniprot
4. Align transcript evidence using Minimap2
5. Merge the predictions from all sources
6. Output the predicted genes to ``predict_results/funannotate_predict.gff3``

Step 3: Functionally Annotate Genes
---------------------------------

The next step is to functionally annotate the predicted genes:

.. code-block:: bash

    funannotate2 annotate --gff3 predict_results/funannotate_predict.gff3 --fasta cleaned_genome.fasta \
        -o annotate_results -s "Aspergillus fumigatus" --strain "Af293" \
        --pfam --dbcan --merops --swissprot --busco --busco_db fungi --cpus 16

This command will:

1. Search the predicted proteins against the Pfam database
2. Search the predicted proteins against the dbCAN database
3. Search the predicted proteins against the MEROPS database
4. Search the predicted proteins against the SwissProt database
5. Search the predicted proteins against the BUSCO database
6. Add the functional annotations to the gene models
7. Output the annotated genes to various formats (GFF3, GenBank, FASTA, etc.)

Step 4: Compare with Another Genome (Optional)
-------------------------------------------

If you have multiple genome annotations, you can compare them:

.. code-block:: bash

    funannotate2 compare -i annotate_results other_genome_results -o compare_results -n "Af293" "Other" --cpus 16

This command will:

1. Compare the gene models between the two genomes
2. Generate various comparison reports
3. Output the comparison results to ``compare_results``

Output Files
-----------

The annotation process produces various output files:

1. **Cleaned Genome**:
   - ``cleaned_genome.fasta``: Cleaned genome assembly

2. **Gene Prediction**:
   - ``predict_results/funannotate_predict.gff3``: Predicted genes in GFF3 format
   - ``predict_results/augustus.gff3``: Augustus predictions
   - ``predict_results/genemark.gtf``: GeneMark predictions
   - ``predict_results/miniprot.gff3``: Miniprot alignments
   - ``predict_results/minimap2_transcripts.gff3``: Minimap2 transcript alignments
   - ``predict_results/minimap2_proteins.gff3``: Minimap2 protein alignments

3. **Functional Annotation**:
   - ``annotate_results/Aspergillus_fumigatus_Af293.gff3``: Annotated genes in GFF3 format
   - ``annotate_results/Aspergillus_fumigatus_Af293.gbk``: Annotated genes in GenBank format
   - ``annotate_results/Aspergillus_fumigatus_Af293.proteins.fa``: Predicted proteins in FASTA format
   - ``annotate_results/Aspergillus_fumigatus_Af293.transcripts.fa``: Predicted transcripts in FASTA format
   - ``annotate_results/Aspergillus_fumigatus_Af293.fasta``: Genome assembly in FASTA format
   - ``annotate_results/Aspergillus_fumigatus_Af293.summary.json``: Summary statistics in JSON format

4. **Genome Comparison**:
   - ``compare_results/summary.tsv``: Summary of the comparison
   - ``compare_results/orthologs.tsv``: Orthologous genes between the genomes
   - ``compare_results/unique_genes.tsv``: Unique genes in each genome
   - ``compare_results/stats.json``: Detailed statistics in JSON format

Troubleshooting
-------------

If you encounter any issues during the annotation process, here are some common solutions:

1. **GeneMark-ES fails**:
   - Make sure GeneMark-ES is installed correctly
   - Check that the genome assembly is not too fragmented
   - Try using a different GeneMark mode (e.g., ET instead of ES)

2. **Augustus fails**:
   - Make sure Augustus is installed correctly
   - Check that the species model exists
   - Try using a different species model

3. **Miniprot/Minimap2 fails**:
   - Make sure Miniprot/Minimap2 is installed correctly
   - Check that the protein/transcript evidence is in the correct format
   - Try using different alignment parameters

4. **Functional annotation fails**:
   - Make sure the required databases are installed correctly
   - Check that the predicted proteins are in the correct format
   - Try using different search parameters

For more help, see the :doc:`faq` or open an issue on the `GitHub repository <https://github.com/nextgenusfs/funannotate2/issues>`_.
