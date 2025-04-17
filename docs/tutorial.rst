Tutorial
========

This tutorial will guide you through the process of annotating a fungal genome using Funannotate2. We will use a publically available genome assembly for Aspergillus nidulans.

Prerequisites
------------

Before starting this tutorial, make sure you have:

1. Installed Funannotate2 and its dependencies and databases (see :doc:`installation`)


Step 1: Fetch genome assembly
--------------------------------

Lets fetch the model organism Aspergillus nidulans FGSCA4 from NCBI:

.. code-block:: bash

    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/011/425/GCF_000011425.1_ASM1142v1/GCF_000011425.1_ASM1142v1_genomic.fna.gz
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/011/425/GCF_000011425.1_ASM1142v1/GCF_000011425.1_ASM1142v1_genomic.gff.gz


Next we can clean up the NCBI GFF3 format:

.. code-block:: bash

   $ gfftk sanitize -f GCF_000011425.1_ASM1142v1_genomic.fna.gz -g GCF_000011425.1_ASM1142v1_genomic.gff.gz -o FGSCA4.gff3


And then we can have a look at the annotation stats in this heavily curated model organism:

.. code-block:: bash

   $ gfftk stats -f GCF_000011425.1_ASM1142v1_genomic.fna.gz -i FGSCA4.gff3
   {
   "genes": 10518,
   "common_name": 112,
   "mRNA": 10455,
   "tRNA": 0,
   "ncRNA": 0,
   "rRNA": 0,
   "repeat_region": 0,
   "misc_feature": 0,
   "avg_gene_length": 1764.49,
   "transcript-level": {
      "CDS_transcripts": 10455,
      "CDS_five_utr": 0,
      "CDS_three_utr": 0,
      "CDS_no_utr": 10455,
      "CDS_five_three_utr": 0,
      "CDS_complete": 10443,
      "CDS_no-start": 3,
      "CDS_no-stop": 7,
      "CDS_no-start_no-stop": 2,
      "total_exons": 35133,
      "total_cds_exons": 34802,
      "average_number_transcripts_per_gene": 0.99,
      "multiple_exon_transcript": 9084,
      "single_exon_transcript": 1371,
      "average_number_cds_exons": 3.33,
      "avg_exon_length": 462.21,
      "median_number_exons": 3,
      "max_number_exons": 26,
      "avg_protein_length": 485.62,
      "avg_transcript_length": 1553.22,
      "functional": {
         "go_terms": 0,
         "interproscan": 8169,
         "eggnog": 0,
         "pfam": 0,
         "cazyme": 0,
         "merops": 0,
         "busco": 0,
         "secretion": 0
      }
   }

Step 2: Train Ab Initio Prediction Tools
-------------------

The next step is to train the ab initio prediction tools, to do that we'll use `funannotate2 train`.  I'm going to use the `anid_f2` as the output directory.

.. code-block:: bash

   $ funannotate2 train -f GCF_000011425.1_ASM1142v1_genomic.fna.gz \
      -g FGSCA4.gff3 -s "Aspergillus nidulans" --strain FGSCA4 --cpus 8 \
      -o anid_f2


.. admonition:: See the full output of this command
    :class: dropdown

      .. code-block:: bash

         [Apr 12 10:12 PM] Python v3.9.19; funannotate2 v25.4.12; gfftk v25.4.12; buscolite v25.4.3
         [Apr 12 10:12 PM] Loading genome assembly and running QC checks
         [Apr 12 10:12 PM] Genome stats:
         {
            "n_contigs": 8,
            "size": 29828291,
            "n50": 3704807,
            "n90": 3145033,
            "l50": 4,
            "l90": 7,
            "avg_length": 3728536
         }
         [Apr 12 10:12 PM] Getting taxonomy information
         {
            "superkingdom": "Eukaryota",
            "kingdom": "Fungi",
            "phylum": "Ascomycota",
            "class": "Eurotiomycetes",
            "order": "Eurotiales",
            "family": "Aspergillaceae",
            "genus": "Aspergillus",
            "species": "Aspergillus nidulans"
         }
         [Apr 12 10:12 PM] Choosing best augustus species based on taxonomy: anidulans
         [Apr 12 10:12 PM] Choosing best busco species based on taxonomy: eurotiales
         [Apr 12 10:12 PM] Running buscolite to generate training set
         [Apr 12 10:12 PM] eurotiales_odb10 lineage contains 4191 BUSCO models
         [Apr 12 10:12 PM] Prefiltering predictions using miniprot of ancestral sequences
         [Apr 12 10:12 PM] Found 725 complete models from miniprot, now launching 3380 augustus/pyhmmer [species=anidulans] jobs for 3379 BUSCO models
         [Apr 12 10:41 PM] Found 3989 BUSCOs in first pass, trying harder to find remaining 202
         [Apr 12 10:41 PM] Found 52 from miniprot, now launching 147 augustus/pyhmmer jobs for 135 BUSCO models
         [Apr 12 10:43 PM] Analysis complete:
            single-copy=4077
            fragmented=0
            duplicated=0
            total=4077
         [Apr 12 10:43 PM] Training set [/Users/jon/software/funannotate2/local_tests/anidulans/train_misc/busco_training_set.gff3] loaded with 4077 gene models
         [Apr 12 10:44 PM] 3,696 of 4,077 models pass training parameters
         [Apr 12 10:44 PM] 3696 gene models selected for training, now splitting into test [n=200] and train [n=3496]
         [Apr 12 10:44 PM] Training augustus using training set
         [Apr 12 10:45 PM] Initial training completed in 00:01:55s
         {
            "tool": "augustus",
            "model": "2729fffa-bec0-45a2-a0fe-b64c0d6ea542",
            "n_test_genes": 200,
            "ref_genes_found": 199,
            "ref_genes_missed": 1,
            "extra_query_genes": 101,
            "average_aed": 0.07467057536626677,
            "nucleotide_sensitivity": 0.9220365983327615,
            "nucleotide_precision": 0.9506290384745041,
            "exon_sensitivity": 0.7030456852791879,
            "exon_precision": 0.7353456611070821,
            "gene_sensitivity": 0.99,
            "gene_precision": 0.495
         }
         [Apr 12 10:45 PM] Training snap using training set
         [Apr 12 10:46 PM] Initial training completed in 00:00:10s
         {
            "tool": "snap",
            "model": "snap-trained.hmm",
            "n_test_genes": 200,
            "ref_genes_found": 200,
            "ref_genes_missed": 0,
            "extra_query_genes": 200,
            "average_aed": 0.11985835682750766,
            "nucleotide_sensitivity": 0.8578286982555101,
            "nucleotide_precision": 0.9623470985417217,
            "exon_sensitivity": 0.5644329896907216,
            "exon_precision": 0.6013132056946491,
            "gene_sensitivity": 1.0,
            "gene_precision": 0.23954372623574144
         }
         [Apr 12 10:46 PM] Training glimmerHMM using training set
         [Apr 12 11:14 PM] Initial training completed in 00:20:17 and parameter optimization completed in 00:07:47s
         {
            "tool": "glimmerhmm",
            "model": "train",
            "n_test_genes": 200,
            "ref_genes_found": 191,
            "ref_genes_missed": 9,
            "extra_query_genes": 90,
            "average_aed": 0.09936167211746938,
            "nucleotide_sensitivity": 0.8940046590916744,
            "nucleotide_precision": 0.9345785751153856,
            "exon_sensitivity": 0.5783783783783784,
            "exon_precision": 0.61981981981982,
            "gene_sensitivity": 0.8846153846153846,
            "gene_precision": 0.4339622641509434
         }
         [Apr 12 11:14 PM] Training GeneMark-ES using self-training
         [Apr 13 02:59 AM] Initial training completed in 03:44:55s
         {
            "tool": "genemark",
            "model": "gmhmm.mod",
            "n_test_genes": 200,
            "ref_genes_found": 200,
            "ref_genes_missed": 0,
            "extra_query_genes": 183,
            "average_aed": 0.062178024762870994,
            "nucleotide_sensitivity": 0.9213744271525245,
            "nucleotide_precision": 0.9748335923946361,
            "exon_sensitivity": 0.745,
            "exon_precision": 0.7820714285714284,
            "gene_sensitivity": 1.0,
            "gene_precision": 0.3879598662207358
         }
         [Apr 13 02:59 AM] Ab initio training finished: /Users/jon/software/funannotate2/local_tests/anidulans/train_results/Aspergillus_nidulans_FGSCA4.params.json
         [Apr 13 02:59 AM] The params.json file can be passed to funannotate2 predict or installed globally with funannotate2 species
         [Apr 13 02:59 AM] funannotate2.train module finished: peak memory usage=204.64 MiB


Step 3: Predict Genes
-------------------

The next step is to predict genes in the cleaned genome assembly:

.. code-block:: bash

    funannotate2 predict -i anid_f2 --cpus 8


.. seealso::
    :class: dropdown

      .. code-block:: bash



This command will:

1. Run GeneMark-ES to predict genes
2. Run Augustus to predict genes
3. Align protein evidence using Miniprot
4. Align transcript evidence using Minimap2
5. Merge the predictions from all sources into consensus models using GFFtk
6. Output the predicted genes to ``predict_results/funannotate_predict.gff3``

Step 4: Functionally Annotate Genes
---------------------------------

The next step is to functionally annotate the predicted genes:

.. code-block:: bash

    funannotate2 annotate -i anidulans --cpus 8

This command will:

1. Search the predicted proteins against the Pfam database
2. Search the predicted proteins against the dbCAN database
3. Search the predicted proteins against the MEROPS database
4. Search the predicted proteins against the SwissProt database
5. Search the predicted proteins against the BUSCO database
6. Add the functional annotations to the gene models
7. Output the annotated genes to various formats (GFF3, GenBank, FASTA, etc.)


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
