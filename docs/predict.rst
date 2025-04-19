Gene Prediction
==============

The ``funannotate2 predict`` command predicts gene models in a eukaryotic genome. It uses a combination of ab initio gene predictors and evidence-based approaches to generate accurate gene models.

Basic Usage
----------

.. code-block:: bash

    funannotate2 predict -f genome.fasta -o predict_results -p pretrained_species -s "Aspergillus fumigatus"

Required Arguments
----------------

* ``-i, --input-dir``: funannotate2 output directory
* ``-f, --fasta``: Input genome FASTA file (softmasked repeats)
* ``-o, --out``: Output folder name
* ``-p, --params, --pretrained``: Params.json or pretrained species slug. Use ``funannotate2 species`` to see pretrained species
* ``-s, --species``: Species name (e.g., "Aspergillus fumigatus")

Optional Arguments
----------------

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

Prediction Process
----------------

The ``predict`` command performs the following steps:

1. **Prepare the genome**:

   * The genome is analyzed for assembly statistics and softmasked regions
   * If the genome is not softmasked, `pytantan <https://github.com/nextgenusfs/pytantan>`_ is used to quickly softmask repeats
   * The genome is split into contigs for parallel processing

2. **Align evidence** (if provided):

   * **Protein evidence**: Proteins are aligned to the genome using `miniprot <https://github.com/lh3/miniprot>`_
   * **Transcript evidence**: Transcripts are aligned to the genome using `gapmm2 <https://github.com/nextgenusfs/gapmm2>`_
   * Evidence alignments are converted to hints for `augustus <https://bioinf.uni-greifswald.de/augustus/>`_

3. **Run ab initio gene predictors**:

   * **`Augustus <https://bioinf.uni-greifswald.de/augustus/>`_**: Uses species-specific parameters and evidence hints
   * **`GeneMark <http://exon.gatech.edu/GeneMark/>`_**: Uses self-training or species-specific parameters
   * **`SNAP <https://github.com/KorfLab/SNAP>`_**: Uses species-specific parameters
   * **`GlimmerHMM <https://ccb.jhu.edu/software/glimmerhmm/>`_**: Uses species-specific parameters
   * **`tRNAscan-SE <http://lowelab.ucsc.edu/tRNAscan-SE/>`_**: Identifies tRNA genes

   Note that you can run other ab initio predictors and they can be passed the :code:`--other_gff` option. For example, `helixerlite <https://github.com/nextgenusfs/helixerlite>`_ is a recommended add-on that uses machine learning for gene prediction, but it's not included in the default installation due to dependency issues.


4. **Generate consensus gene models**:

   * The `GFFtk <https://github.com/nextgenusfs/gfftk>`_ consensus module is used to integrate all evidence and ab initio predictions
   * Gene models are weighted based on the reliability of each source
   * Overlapping gene models are resolved based on evidence and prediction quality
   * Gene models are filtered based on various criteria (e.g., minimum protein length)

5. **Annotate gene models**:

   * Gene models are assigned unique IDs based on the locus tag and numbering
   * Gene models are sorted by genomic location
   * Gene models are output in GFF3, TBL, and GenBank formats
   * Protein and transcript sequences are extracted from the gene models
   * Summary statistics are generated for the annotation

Evidence-Based Prediction
----------------------

For best results, provide protein and/or transcript evidence. By defualt, the UniProt/SwissProt database is used for protein evidence.

.. code-block:: bash

    funannotate2 predict -f genome.fasta -o predict_results -p pretrained_species -s "Aspergillus fumigatus" \
        -ps uniprot_fungi.fasta -ts rnaseq_transcripts.fasta

Protein evidence should be in FASTA format and can include:

* Proteins from closely related species
* Curated protein databases (e.g., UniProt)
* Proteins from previous annotations

Transcript evidence should be in FASTA format and can include:

* Assembled transcripts from RNA-Seq data
* EST sequences
* cDNA sequences

Using External Gene Models
-----------------------

You can provide external gene models in GFF3 format:

.. code-block:: bash

    funannotate2 predict -f genome.fasta -o f2_output -p pretrained_species -s "Aspergillus fumigatus" -e external_models.gff3

External gene models can be from:

* Previous annotations
* Other gene prediction tools
* Manual annotations

Output Files
----------

The ``predict`` command generates the following output files in the specified output directory:

* **<species>.gff3**: Gene models in GFF3 format
* **<species>.tbl**: Gene models in NCBI TBL format
* **<species>.gbk**: Gene models in GenBank format
* **<species>.proteins.fa**: Protein sequences in FASTA format
* **<species>.transcripts.fa**: Transcript sequences in FASTA format
* **<species>.fasta**: Genome sequence in FASTA format
* **<species>.summary.json**: Summary statistics in JSON format

The ``predict_misc`` directory contains intermediate files and detailed results from each prediction source:

* **augustus/**: Augustus prediction results
* **genemark/**: GeneMark prediction results
* **snap/**: SNAP prediction results
* **glimmerhmm/**: GlimmerHMM prediction results
* **trnascan/**: tRNAscan-SE results
* **proteins/**: Protein evidence alignments
* **transcripts/**: Transcript evidence alignments
* **hints/**: Evidence hints for ab initio predictors
* **softmasked-regions.bed**: Softmasked regions in BED format
* **assembly-gaps.bed**: Assembly gaps in BED format
