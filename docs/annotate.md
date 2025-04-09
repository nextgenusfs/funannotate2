# Functional Annotation

The `funannotate2 annotate` command adds functional annotation to gene models. It can use various sources of annotation, including InterProScan, EggNOG-mapper, antiSMASH, and more.

## Basic Usage

```bash
funannotate2 annotate -i /path/to/funannotate2_predict_output -o /path/to/output_dir
```

## Required Arguments

- `-i, --input-dir`: Path to funannotate2 predict output directory
- `-f, --fasta`: Genome in FASTA format (required if not using --input-dir)
- `-t, --tbl`: Genome annotation in TBL format (required if not using --input-dir and not using --gff3)
- `-g, --gff3`: Genome annotation in GFF3 format (required if not using --input-dir and not using --tbl)
- `-o, --out`: Output folder name (required if not using --input-dir)

## Optional Arguments

- `-a, --annotations`: Annotations files, 3 column TSV [transcript-id, feature, data]
- `-s, --species`: Species name, use quotes for binomial, e.g. "Aspergillus fumigatus"
- `-st, --strain`: Strain/isolate name
- `--cpus`: Number of CPUs to use (default: 2)
- `--tmpdir`: Volume to write tmp files (default: /tmp)
- `--curated-names`: Path to custom file with gene-specific annotations (tab-delimited: gene_id\tannotation_type\tannotation_value)

## Using Custom Curated Gene Names and Products

Funannotate2 includes a database of curated gene names and products that are used to ensure consistent and accurate annotation. However, you may want to use your own curated gene names and products for specific organisms or projects.

### Creating a Custom Annotations File

Create a tab-delimited text file with gene IDs in the first column, annotation types in the second column, and annotation values in the third column:

```
# Custom annotations for specific genes/transcripts
gene123	name	ACT1
gene123	product	Actin
gene456	name	CDC42
gene456	product	Cell division control protein 42
gene789	go_term	GO:0005524
```

A template file is available at `funannotate2/data/custom_annotations.template.txt`.

### Using the Custom Annotations File

Use the `--curated-names` option to specify your custom file:

```bash
funannotate2 annotate -i /path/to/funannotate2_predict_output -o /path/to/output_dir --curated-names /path/to/custom_annotations.txt
```

### How Custom Annotations are Used

1. Funannotate2 will load your custom annotations file
2. For each gene ID in your custom file, the specified annotations will be applied
3. For single-value annotations like `name` and `product`, custom values replace existing ones
4. For multi-value annotations like `go_term` and `ec_number`, custom values are added to existing ones
5. Custom annotations for `name` and `product` bypass the cleaning rules
6. This allows you to have precise control over important annotations while preserving other data

### Supported Annotation Types

- `name`: Gene name (e.g., ACT1, CDC42)
- `product`: Product description (e.g., Actin, Cell division control protein 42)
- `note`: Additional information about the gene
- `go_term`: Gene Ontology term (e.g., GO:0005524)
- `ec_number`: Enzyme Commission number (e.g., 3.6.4.13)
- `db_xref`: Database cross-reference (e.g., UniProtKB:P60010)

### Benefits of Using Custom Annotations

- Precise control over annotations for specific genes
- Override automated annotation for important genes
- Add specialized annotations not available from automated sources
- Ensure consistent annotation across multiple projects
- Maintain control over the annotation of important genes
