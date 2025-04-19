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


And then we can have a look at the annotation stats in this heavily curated model organism, and you can see there are 10518 genes, 10455 of them are protein coding or (mRNA) features.

.. code-block:: bash

   $ gfftk stats -f GCF_000011425.1_ASM1142v1_genomic.fna.gz -i FGSCA4.gff3



.. admonition:: See the full output of this command
    :class: dropdown

      .. code-block:: bash

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

The next step is to train the ab initio prediction tools, to do that we'll use ``funannotate2 train``.  I'm going to use the ``anid_f2`` as the output directory.

.. code-block:: bash

   $ funannotate2 train -f GCF_000011425.1_ASM1142v1_genomic.fna.gz \
      -s "Aspergillus nidulans" --strain FGSCA4 --cpus 8 -o anid_f2


.. admonition:: See the full output of this command
    :class: dropdown

      .. code-block:: bash

         [Apr 16 07:53 PM] Python v3.9.19; funannotate2 v25.4.15; gfftk v25.4.12; buscolite v25.4.3
         [Apr 16 07:53 PM] Loading genome assembly and running QC checks
         [Apr 16 07:53 PM] Genome stats:
         {
            "n_contigs": 8,
            "size": 29828291,
            "n50": 3704807,
            "n90": 3145033,
            "l50": 4,
            "l90": 7,
            "avg_length": 3728536
         }
         [Apr 16 07:53 PM] Getting taxonomy information
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
         [Apr 16 07:53 PM] Choosing best augustus species based on taxonomy: aspergillus_nidulans
         [Apr 16 07:53 PM] Choosing best busco species based on taxonomy: eurotiales
         [Apr 16 07:53 PM] Running buscolite to generate training set
         [Apr 16 07:53 PM] eurotiales_odb10 lineage contains 4191 BUSCO models
         [Apr 16 07:53 PM] Prefiltering predictions using miniprot of ancestral sequences
         [Apr 16 07:53 PM] Found 725 complete models from miniprot, now launching 3380 augustus/pyhmmer [species=aspergillus_nidulans] jobs for 3379 BUSCO models
         [Apr 16 08:20 PM] Found 3984 BUSCOs in first pass, trying harder to find remaining 207
         [Apr 16 08:20 PM] Found 50 from miniprot, now launching 154 augustus/pyhmmer jobs for 142 BUSCO models
         [Apr 16 08:22 PM] Analysis complete:
            single-copy=4071
            fragmented=0
            duplicated=0
            total=4071
         [Apr 16 08:22 PM] Training set [/Users/jon/software/funannotate2/local_tests/anid_f2/train_misc/busco_training_set.gff3] loaded with 4071 gene models
         [Apr 16 08:22 PM] 3,645 of 4,071 models pass training parameters
         [Apr 16 08:22 PM] 3645 gene models selected for training, now splitting into test [n=200] and train [n=3445]
         [Apr 16 08:22 PM] Training augustus using training set
         [Apr 16 08:24 PM] Initial training completed in 00:01:50s
         {
            "tool": "augustus",
            "model": "3edf7f79-7a02-4685-b626-90dbadd56e97",
            "n_test_genes": 200,
            "ref_genes_found": 196,
            "ref_genes_missed": 4,
            "extra_query_genes": 104,
            "average_aed": 0.0730729493547641,
            "nucleotide_sensitivity": 0.919880412361861,
            "nucleotide_precision": 0.954487814565406,
            "exon_sensitivity": 0.7239583333333334,
            "exon_precision": 0.77092539983165,
            "gene_sensitivity": 0.9622641509433962,
            "gene_precision": 0.49514563106796117
         }
         [Apr 16 08:24 PM] Training snap using training set
         [Apr 16 08:25 PM] Initial training completed in 00:00:10s
         {
            "tool": "snap",
            "model": "snap-trained.hmm",
            "n_test_genes": 200,
            "ref_genes_found": 200,
            "ref_genes_missed": 0,
            "extra_query_genes": 184,
            "average_aed": 0.12570097490993856,
            "nucleotide_sensitivity": 0.8592611345633322,
            "nucleotide_precision": 0.9431127619612985,
            "exon_sensitivity": 0.5634517766497462,
            "exon_precision": 0.6074208363548468,
            "gene_sensitivity": 1.0,
            "gene_precision": 0.2550607287449393
         }
         [Apr 16 08:25 PM] Training glimmerHMM using training set
         [Apr 16 08:53 PM] Initial training completed in 00:21:03 and parameter optimization completed in 00:07:48s
         {
            "tool": "glimmerhmm",
            "model": "train",
            "n_test_genes": 200,
            "ref_genes_found": 199,
            "ref_genes_missed": 1,
            "extra_query_genes": 140,
            "average_aed": 0.10988379675849991,
            "nucleotide_sensitivity": 0.8847082643612956,
            "nucleotide_precision": 0.9255094638711031,
            "exon_sensitivity": 0.6395939086294417,
            "exon_precision": 0.6776669889614054,
            "gene_sensitivity": 0.9873417721518988,
            "gene_precision": 0.3577981651376147
         }
         [Apr 16 08:53 PM] Training GeneMark-ES using self-training
         [Apr 16 11:31 PM] Initial training completed in 02:37:43s
         {
            "tool": "genemark",
            "model": "gmhmm.mod",
            "n_test_genes": 200,
            "ref_genes_found": 200,
            "ref_genes_missed": 0,
            "extra_query_genes": 180,
            "average_aed": 0.0625898722361626,
            "nucleotide_sensitivity": 0.9294058289674456,
            "nucleotide_precision": 0.9638773813305964,
            "exon_sensitivity": 0.76,
            "exon_precision": 0.8057947330447329,
            "gene_sensitivity": 1.0,
            "gene_precision": 0.3939393939393939
         }
         [Apr 16 11:31 PM] Ab initio training finished: /Users/jon/software/funannotate2/local_tests/anid_f2/train_results/Aspergillus_nidulans_FGSCA4.params.json
         [Apr 16 11:31 PM] The params.json file can be passed to funannotate2 predict or installed globally with funannotate2 species
         [Apr 16 11:31 PM] funannotate2.train module finished: peak memory usage=203.82 MiB


Step 3: Predict Genes
-------------------

The next step is to predict genes using the training sets we just generated.  Here we will just use the defaults for evidence mapping, which is to align the SwissProt/UniProt curated database. If you had high quality curated transcript data you could pass that to ``--transcripts`` option.

.. code-block:: bash

   $ funannotate2 predict -i anid_f2 --cpus 8 --strain FGSCA4


.. admonition:: See the full output of this command
    :class: dropdown

      .. code-block:: bash

         [Apr 17 11:16 AM] Python v3.9.19; funannotate2 v25.4.15; gfftk v25.4.13; buscolite v25.4.3
         [Apr 17 11:16 AM] Parsed data from --input-dir anid_f2
            --fasta /Users/jon/software/funannotate2/local_tests/anid_f2/train_results/GCF_000011425.1_ASM1142v1_genomic.fna.gz
            --species "Aspergillus nidulans"
            --params /Users/jon/software/funannotate2/local_tests/anid_f2/train_results/Aspergillus_nidulans_FGSCA4.params.json
            --out anid_f2
         [Apr 17 11:16 AM] Loaded training params for Aspergillus_nidulans_FGSCA4: ['augustus', 'glimmerhmm', 'snap', 'genemark']
         [Apr 17 11:16 AM] temporary files located in: /tmp/predict_d61fc66a-78ea-480a-a385-eb7e2c71ab88
         [Apr 17 11:16 AM] Loading genome assembly, running QC checks, searching for mitochondrial contigs, calculating softmasked regions and assembly gaps
         [Apr 17 11:17 AM] Genome stats:
         {
            "n_contigs": 8,
            "size": 29828291,
            "softmasked": "5.10%",
            "gaps": "0.03%",
            "n50": 3704807,
            "n90": 3145033,
            "l50": 4,
            "l90": 7,
            "avg_length": 3728536
         }
         [Apr 17 11:17 AM] Parsed 484021 [out of 572214] proteins to align for evidence, derived from:
         /Users/jon/software/funannotate_db/uniprot_sprot.fasta
         [Apr 17 11:17 AM] Aligning protein evidence to the genome assembly with miniprot
         [Apr 17 11:37 AM] Generated 100013 alignments: 824 were valid gene models
         [Apr 17 11:37 AM] Parsing alignments and generating hintsfile for augustus
         [Apr 17 11:37 AM] Running ab initio gene predictions using 8 cpus
         [Apr 17 01:08 PM] Ab initio predictions finished:
         {
            "augustus": 7048,
            "augustus-hiq": 1480,
            "glimmerhmm": 9195,
            "snap": 10138,
            "genemark": 10379
         }
         [Apr 17 01:08 PM] Measuring assembly completeness with buscolite for all ab initio predictions
         [Apr 17 01:11 PM] ab initio models scoring by algorithm:
         {
            "augustus": {
               "busco": 1.0065963060686016,
               "train": {
                  "tool": "augustus",
                  "model": "3edf7f79-7a02-4685-b626-90dbadd56e97",
                  "n_test_genes": 200,
                  "ref_genes_found": 196,
                  "ref_genes_missed": 4,
                  "extra_query_genes": 104,
                  "average_aed": 0.0730729493547641,
                  "nucleotide_sensitivity": 0.919880412361861,
                  "nucleotide_precision": 0.954487814565406,
                  "exon_sensitivity": 0.7239583333333334,
                  "exon_precision": 0.77092539983165,
                  "gene_sensitivity": 0.9622641509433962,
                  "gene_precision": 0.49514563106796117
               }
            },
            "genemark": {
               "busco": 0.9947229551451188,
               "train": {
                  "tool": "genemark",
                  "model": "gmhmm.mod",
                  "n_test_genes": 200,
                  "ref_genes_found": 200,
                  "ref_genes_missed": 0,
                  "extra_query_genes": 180,
                  "average_aed": 0.0625898722361626,
                  "nucleotide_sensitivity": 0.9294058289674456,
                  "nucleotide_precision": 0.9638773813305964,
                  "exon_sensitivity": 0.76,
                  "exon_precision": 0.8057947330447329,
                  "gene_sensitivity": 1.0,
                  "gene_precision": 0.3939393939393939
               }
            },
            "glimmerhmm": {
               "busco": 0.9736147757255936,
               "train": {
                  "tool": "glimmerhmm",
                  "model": "train",
                  "n_test_genes": 200,
                  "ref_genes_found": 199,
                  "ref_genes_missed": 1,
                  "extra_query_genes": 140,
                  "average_aed": 0.10988379675849991,
                  "nucleotide_sensitivity": 0.8847082643612956,
                  "nucleotide_precision": 0.9255094638711031,
                  "exon_sensitivity": 0.6395939086294417,
                  "exon_precision": 0.6776669889614054,
                  "gene_sensitivity": 0.9873417721518988,
                  "gene_precision": 0.3577981651376147
               }
            },
            "miniprot-gene": {
               "busco": 0.158311345646438
            },
            "snap": {
               "busco": 0.9472295514511874,
               "train": {
                  "tool": "snap",
                  "model": "snap-trained.hmm",
                  "n_test_genes": 200,
                  "ref_genes_found": 200,
                  "ref_genes_missed": 0,
                  "extra_query_genes": 184,
                  "average_aed": 0.12570097490993856,
                  "nucleotide_sensitivity": 0.8592611345633322,
                  "nucleotide_precision": 0.9431127619612985,
                  "exon_sensitivity": 0.5634517766497462,
                  "exon_precision": 0.6074208363548468,
                  "gene_sensitivity": 1.0,
                  "gene_precision": 0.2550607287449393
               }
            }
         }
         [Apr 17 01:11 PM] Calculated ab initio weights from data: ['augustus:2', 'genemark:3', 'glimmerhmm:1', 'miniprot-gene:1', 'snap:1', 'augustus-hiq:4']
         [Apr 17 01:11 PM] GFFtk consensus will generate the best gene model at each locus
         [Apr 17 01:11 PM] Parsing GFF3 files and clustering data into strand specific loci
         [Apr 17 01:13 PM] Merging gene predictions from 6 source files
            {'predictions': {'augustus-hiq', 'snap', 'augustus', 'genemark', 'miniprot-gene', 'glimmerhmm'}, 'evidence': {'miniprot'}}
         [Apr 17 01:13 PM] Parsed 39064 gene models into 8708 loci. Dropped 81 genes models that were pseudo [labled as such or internal stop codons]
            {'augustus': 81}
         [Apr 17 01:13 PM] Using these filtered loci, the calculated gene model source weights to use as tiebreakers:
            {"augustus-hiq": 4, "genemark": 3, "augustus": 2, "snap": 1, "glimmerhmm": 1, "miniprot-gene": 1}
         [Apr 17 01:13 PM] Processing 8708 loci using 8 processes
         [Apr 17 01:13 PM] Setting minimum gene model score to 7
         [Apr 17 01:13 PM] Loaded repeats representing 5.10% of the genome and filtering out loci that are > 90% overlap with repeats
         [Apr 17 01:13 PM] 60 gene models were dropped due to repeat overlap
         [Apr 17 01:13 PM] 10016 consensus gene models derived from these sources:
            [["genemark", 8235], ["augustus-hiq", 1400], ["augustus", 164], ["glimmerhmm", 154], ["snap", 122], ["miniprot-gene", 1]]
         [Apr 17 01:13 PM] GFFtk consensus is finished: /Users/jon/software/funannotate2/local_tests/anid_f2/predict_misc/consensus.predictions.gff3
         [Apr 17 01:13 PM] Predicting tRNA genes
         [Apr 17 01:13 PM] Merging all gene models, sorting, and renaming using locus_tag=FUN2_
         [Apr 17 01:14 PM] Converting to GenBank format
         [Apr 17 01:14 PM] Annotation statistics:
         {
            "genes": 10196,
            "common_name": 0,
            "mRNA": 10016,
            "tRNA": 180,
            "ncRNA": 0,
            "rRNA": 0,
            "repeat_region": 0,
            "misc_feature": 0,
            "avg_gene_length": 1658.12,
            "transcript-level": {
               "CDS_transcripts": 10016,
               "CDS_five_utr": 0,
               "CDS_three_utr": 0,
               "CDS_no_utr": 10016,
               "CDS_five_three_utr": 0,
               "CDS_complete": 10012,
               "CDS_no-start": 2,
               "CDS_no-stop": 2,
               "CDS_no-start_no-stop": 0,
               "total_exons": 30954,
               "total_cds_exons": 30954,
               "average_number_transcripts_per_gene": 1.0,
               "multiple_exon_transcript": 7892,
               "single_exon_transcript": 2124,
               "average_number_cds_exons": 3.09,
               "avg_exon_length": 490.66,
               "median_number_exons": 3.0,
               "max_number_exons": 27,
               "avg_protein_length": 505.49,
               "avg_transcript_length": 1516.37,
               "functional": {
                  "go_terms": 0,
                  "interproscan": 0,
                  "eggnog": 0,
                  "pfam": 0,
                  "cazyme": 0,
                  "merops": 0,
                  "busco": 0,
                  "secretion": 0
               }
            }
         }
         [Apr 17 01:15 PM] Measuring assembly completeness with buscolite
         [Apr 17 01:15 PM] Assembly completeness:
            complete=750 [98.94%]
            single-copy=693 [91.42%]
            fragmented=0 [0.00%]
            duplicated=57 [7.52%]
            missing=8 [1.06%]
            total=758 [100.00%]
         [Apr 17 01:15 PM] funannotate2.predict module finished: peak memory usage=1.11 GiB


In this case we see that the tool predicted 10196 genes (10016 of which are protein coding genes), note this very close the the public gold standard annotation in terms of numbers of gene models. You can see at the end of predict that the tool also ran BUSCO to measure the completeness of the assembly, which is 98.94% complete.

Step 4: Functionally Annotate Genes
---------------------------------

The next step is to functionally annotate the predicted gene models.  The core ``annotate`` module in ``funannotate2`` will do a few functional annotation steps, however, it does not natively try to support everything as this becomes daunting. Rather users can provide their own functional annotation to the script via the ``-a,--annotations`` command line argument.

.. code-block:: bash

   $ funannotate2 annotate -i anid_f2/ --cpus 8


.. admonition:: See the full output of this command
    :class: dropdown

      .. code-block:: bash

         [Apr 17 01:43 PM] Python v3.9.19; funannotate2 v25.4.15; gfftk v25.4.13; buscolite v25.4.3
         [Apr 17 01:43 PM] Parsed input files from --input-dir anid_f2/
         --fasta /Users/jon/software/funannotate2/local_tests/anid_f2/predict_results/Aspergillus_nidulans.fasta
         --tbl /Users/jon/software/funannotate2/local_tests/anid_f2/predict_results/Aspergillus_nidulans.tbl
         --gff3 /Users/jon/software/funannotate2/local_tests/anid_f2/predict_results/Aspergillus_nidulans.gff3
         --out anid_f2/
         [Apr 17 01:43 PM] temporary files located in: /tmp/annotate_c76039d4-7e25-4511-931a-48a8a2b6c7bc
         [Apr 17 01:43 PM] Parsed genome stats:
         [Apr 17 01:43 PM]
         {
            "genes": 10196,
            "common_name": 0,
            "mRNA": 10016,
            "tRNA": 180,
            "ncRNA": 0,
            "rRNA": 0,
            "repeat_region": 0,
            "misc_feature": 0,
            "avg_gene_length": 1658.12,
            "transcript-level": {
               "CDS_transcripts": 10016,
               "CDS_five_utr": 0,
               "CDS_three_utr": 0,
               "CDS_no_utr": 10016,
               "CDS_five_three_utr": 0,
               "CDS_complete": 10012,
               "CDS_no-start": 2,
               "CDS_no-stop": 2,
               "CDS_no-start_no-stop": 0,
               "total_exons": 30954,
               "total_cds_exons": 30954,
               "average_number_transcripts_per_gene": 1.0,
               "multiple_exon_transcript": 7892,
               "single_exon_transcript": 2124,
               "average_number_cds_exons": 3.09,
               "avg_exon_length": 490.66,
               "median_number_exons": 3.0,
               "max_number_exons": 27,
               "avg_protein_length": 505.49,
               "avg_transcript_length": 1516.37,
               "functional": {
                  "go_terms": 0,
                  "interproscan": 0,
                  "eggnog": 0,
                  "pfam": 0,
                  "cazyme": 0,
                  "merops": 0,
                  "busco": 0,
                  "secretion": 0
               }
            }
         }
         [Apr 17 01:43 PM] Annotating proteome with pyhmmer against the Pfam-A database
         [Apr 17 01:46 PM] Pfam-A search resulted in 16464 hits and finished in 137.72 seconds
         [Apr 17 01:46 PM] Annotating proteome with pyhmmer against the dbCAN (CAZyme) database
         [Apr 17 01:46 PM] dbCAN search resulted in 545 hits and finished in 16.15 seconds
         [Apr 17 01:46 PM] Annotating proteome with diamond against the UniProtKB/Swiss-Prot database
         [Apr 17 01:47 PM] UniProtKB/Swiss-Prot search resulted in 1838 hits and finished in 71.78 seconds
         [Apr 17 01:47 PM] Annotating proteome with diamond against the MEROPS protease database
         [Apr 17 01:47 PM] MEROPS search resulted in 313 hits and finished in 1.9 seconds
         [Apr 17 01:47 PM] BUSCOlite [conserved ortholog] search using eurotiales models
         [Apr 17 01:49 PM] BUSCOlite search resulted in 4194 hits and finished in 98.32 seconds
         [Apr 17 01:49 PM] Found functional annotation for 8594 gene models
         [Apr 17 01:49 PM] Annotation sources: {'db_xref': 9903, 'note': 5940, 'name': 925, 'product': 925}
         [Apr 17 01:49 PM] Cleaning gene names and product descriptions using curated database
         [Apr 17 01:49 PM] Found 38 new valid gene names/products that could be added to the curated database
         [Apr 17 01:49 PM] See /Users/jon/software/funannotate2/local_tests/anid_f2/annotate_results/Gene2Products.new-valid.txt for details
         [Apr 17 01:49 PM] Converting to GenBank format
         [Apr 17 01:49 PM] Writing rest of the output annotation files
         [Apr 17 01:49 PM] Annotation Summary:
         [Apr 17 01:49 PM]
         {
            "genes": 10196,
            "common_name": 925,
            "mRNA": 10016,
            "tRNA": 180,
            "ncRNA": 0,
            "rRNA": 0,
            "repeat_region": 0,
            "misc_feature": 0,
            "avg_gene_length": 1658.12,
            "transcript-level": {
               "CDS_transcripts": 10016,
               "CDS_five_utr": 0,
               "CDS_three_utr": 0,
               "CDS_no_utr": 10016,
               "CDS_five_three_utr": 0,
               "CDS_complete": 10012,
               "CDS_no-start": 2,
               "CDS_no-stop": 2,
               "CDS_no-start_no-stop": 0,
               "total_exons": 30954,
               "total_cds_exons": 30954,
               "average_number_transcripts_per_gene": 1.0,
               "multiple_exon_transcript": 7892,
               "single_exon_transcript": 2124,
               "average_number_cds_exons": 3.09,
               "avg_exon_length": 490.66,
               "median_number_exons": 3.0,
               "max_number_exons": 27,
               "avg_protein_length": 505.49,
               "avg_transcript_length": 1516.37,
               "functional": {
                  "go_terms": 0,
                  "interproscan": 0,
                  "eggnog": 0,
                  "pfam": 8065,
                  "cazyme": 520,
                  "merops": 313,
                  "busco": 4194,
                  "secretion": 0
               }
            }
         }
         [Apr 17 01:49 PM] funannotate2.annotate module finished: peak memory usage=183.47 MiB

The annotate step was able to add functional annotation to 8594 of the 10196 gene models from PFAM, CAZymes, MEROPS, SwissProt/UniProt, and BUSCO.

Output Files
-----------

The annotation process produces various output in the output directory (``anid_f2``). Here are some of the key files:

1. **Ab initio training parameters**:

   - ``train_results/Aspergillus_nidulans_FGSCA4.params.json``: Training parameters for ab initio gene prediction tools, can be used for future predictions or permanently install with `funannotate2 species`


2. **Gene Prediction**:

   * ``predict_results/Aspergillus_nidulans_FGSCA4.fasta``: Genome assembly in FASTA format

   * ``predict_results/Aspergillus_nidulans_FGSCA4.gff3``: Predicted genes in GFF3 format

   * ``predict_results/Aspergillus_nidulans_FGSCA4.gbk``: Predicted genes in GenBank flat-file format

   * ``predict_results/Aspergillus_nidulans_FGSCA4.tbl``: Predicted genes in GenBank TBL format


3. **Functional Annotation**:

   * ``annotate_results/Aspergillus_nidulans_FGSCA4.fasta``: Genome assembly in FASTA format

   * ``annotate_results/Aspergillus_nidulans_FGSCA4.gff3``: Predicted genes in GFF3 format

   * ``annotate_results/Aspergillus_nidulans_FGSCA4.gbk``: Predicted genes in GenBank flat-file format

   * ``annotate_results/Aspergillus_nidulans_FGSCA4.tbl``: Predicted genes in GenBank TBL format

   * ``annotate_results/Aspergillus_nidulans_FGSCA4.proteins.fa``: Predicted proteins in FASTA format

   * ``annotate_results/Aspergillus_nidulans_FGSCA4.transcripts.fa``: Predicted transcripts in FASTA format

   * ``annotate_results/Aspergillus_nidulans_FGSCA4.summary.json``: Summary statistics in JSON format



For more help, see the :doc:`faq` or open an issue on the `GitHub repository <https://github.com/nextgenusfs/funannotate2/issues>`_.
