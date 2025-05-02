import os
import sys
import uuid

from gapmm2.align import aligner as transcript_aligner
from gfftk.gff import dict2gff3, dict2gff3alignments, gff2dict

from .config import env
from .utilities import checkfile, execute, merge_coordinates, runSubprocess


def split_evidence_and_genes(gff, fasta, evidence, genes, gff_format="alignment"):
    """
    Parse, validate, and split a GFF3 file into gene models and alignments.

    This function processes a GFF3 file to extract gene models and alignments, filters gene models based on specific criteria, and writes the results to separate files. It ensures that the evidence is in an EVM-compatible format and handles different output formats based on the input GFF format.

    Args:
        gff (str): Path to the input GFF3 file.
        fasta (str): Path to the corresponding FASTA file.
        evidence (str): Output file path for the evidence alignments.
        genes (str): Output file path for the gene models.
        gff_format (str, optional): Format of the input GFF file. Defaults to "alignment".

    Returns:
        tuple: A tuple containing:
            - int: The count of total alignments.
            - int: The count of unique gene models.
    """
    # function to take a GFF3 file, parse, validate, and split
    aligns = gff2dict(gff, fasta, gff_format=gff_format)
    g = {}
    gseen = set()
    for k, v in aligns.items():
        if "mRNA" in v["type"]:
            if (
                False in v["partialStart"]
                and False in v["partialStop"]
                and v["pseudo"] is False
            ):
                # we don't want multiple of the same gene model
                unique_str = f"{v['contig']}_{v['strand']}_{v['location'][0]}_{v['location'][1]}_{len(v['CDS'][0])}"
                if unique_str not in gseen:
                    g[k] = v
                    gseen.add(unique_str)
    # now we can write these to files
    if gff_format == "miniprot":
        new_source = "miniprot-gene"
    else:
        new_source = "gapmm2-gene"
    dict2gff3(g, output=genes, source=new_source)
    # evidence needs to be in EVM compatible format
    # if gff_format == 'alignment':
    #    dict2gff3alignments(aligns, output=evidence, alignments="transcript")
    # if proteins, also need to write alignment file?
    if gff_format == "miniprot":
        dict2gff3alignments(aligns, output=evidence, alignments="protein")
    return len(aligns), len(g)


def align_transcripts(
    fasta,
    transcripts,
    evidence,
    genes,
    cpus=1,
    max_intron=3000,
    log=sys.stderr.write,
    tmpdir="/tmp",
):
    """
    Run spliced transcript alignment using gapmm2 (mm2 + edlib refinement).

    Aligns transcript evidence to the genome assembly with gapmm2, generates an EVM compatible GFF alignment file,
    parses data, and extracts full-length gene models.

    Args:
        fasta (str): Path to the FASTA file of the genome assembly.
        transcripts (str): Path to the transcript evidence file.
        evidence (str): Path to store the aligned evidence in GFF3 format.
        genes (str): Path to store the extracted full-length gene models in GFF3 format.
        cpus (int, optional): Number of CPUs to use for alignment. Defaults to 1.
        max_intron (int, optional): Maximum intron length allowed. Defaults to 3000.
        log (function, optional): Logging function. Defaults to sys.stderr.write.
        tmpdir (str, optional): Temporary directory path. Defaults to "/tmp".
    """
    # function to run spliced transcript alignment using gapmm2 (mm2 + edlib refinement)
    # generate an EVM compatible GFF alignment file, parse data and pull out any full length gene models
    log.info("Aligning transcript evidence to the genome assembly with gapmm2")

    t_stats = transcript_aligner(
        fasta,
        transcripts,
        output=evidence,
        out_fmt="gff3",
        threads=cpus,
        max_intron=max_intron,
    )

    # now we need to parse this evidence GFF3 file and look for full length
    n_aligns, n_genes = split_evidence_and_genes(
        evidence, fasta, evidence, genes, gff_format="alignment"
    )
    log.info(
        f"Generated {t_stats['n']} alignments: {t_stats['refine-left']} required 5' refinement, {t_stats['refine-right']} required 3' refinement, {t_stats['low-mapq']} dropped low score, {n_genes} were valid gene models"
    )


def align_proteins(
    fasta,
    proteins,
    evidence,
    genes,
    cpus=1,
    max_intron=3000,
    log=sys.stderr.write,
    tmpdir="/tmp",
):
    """
    Align protein evidence to the genome assembly using miniprot.

    This function generates evidence alignments and parses full-length models from the alignment results.

    Args:
        fasta (str): Path to the genome assembly FASTA file.
        proteins (str): Path to the protein evidence file.
        evidence (str): Output file path for the evidence alignments.
        genes (str): Output file path for the valid gene models.
        cpus (int, optional): Number of CPUs to use for alignment. Defaults to 1.
        max_intron (int, optional): Maximum intron length for alignment. Defaults to 3000.
        log (function, optional): Logger function for info and error messages. Defaults to sys.stderr.write.
        tmpdir (str, optional): Temporary directory path for storing intermediate files. Defaults to "/tmp".

    Returns:
        None
    """
    # function to align protein evidence to genome with miniprot
    # generate evidence alignments and then parse any full length models
    mini_tmp = os.path.join(tmpdir, f"{uuid.uuid1()}.miniprot.out")
    # miniprot --outn=4 --gff-only -t 8 valid_contigs.fasta /usr/local/share/funannotate/uniprot_sprot.fasta > uniprot.mini.gff3
    log.info("Aligning protein evidence to the genome assembly with miniprot")
    cmd = [
        "miniprot",
        "--outn=4",
        "-G",
        str(max_intron),
        "--gff-only",
        "-t",
        str(cpus),
        os.path.abspath(fasta),
        os.path.abspath(proteins),
    ]
    try:
        runSubprocess(cmd, log, stdout=mini_tmp)
    except Exception as e:
        log.error(f"Error occurred during `runSubprocess`: {e}")
        os.remove(mini_tmp)
    # miniprot uses non-standard gff, so load with gfftk parse
    # now we need to parse this evidence GFF3 file and look for full length
    n_aligns, n_genes = split_evidence_and_genes(
        mini_tmp, fasta, evidence, genes, gff_format="miniprot"
    )
    log.info(f"Generated {n_aligns} alignments: {n_genes} were valid gene models")
    os.remove(mini_tmp)


def align_mito(query, cpus=1, min_cov=0.65, min_qual=50, debug=False):
    """
    Align a query genome against the RefSeq mitochondrial database to identify putative mitochondrial contigs.

    This function uses `minimap2` to perform the alignment and processes the results to calculate coverage.

    Args:
        query (str): Path to the query genome file.
        cpus (int, optional): Number of CPUs to use for alignment. Defaults to 1.
        min_cov (float, optional): Minimum coverage threshold for alignment. Defaults to 0.65.
        min_qual (int, optional): Minimum alignment quality score. Defaults to 50.

    Returns:
        tuple: A tuple containing two dictionaries:
            - The first dictionary contains contigs with sufficient coverage.
            - The second dictionary contains all alignment details.
    """
    # align genome against refseq mitochondrial database to find putative mitochondrial contigs
    d = {}
    f = {}
    mito_db = os.path.join(env["FUNANNOTATE2_DB"], "mito.mmi")
    cmd = [
        "minimap2",
        "-t",
        str(cpus),
        os.path.abspath(mito_db),
        os.path.abspath(query),
    ]
    if checkfile(mito_db):
        for line in execute(cmd):
            line = line.strip()
            data = line.split("\t")
            if int(data[11]) >= min_qual:
                if data[0] not in d:
                    d[data[0]] = {
                        "length": int(data[1]),
                        "alignments": {data[5]: [(int(data[2]), int(data[3]))]},
                    }
                else:
                    if data[5] not in d[data[0]]["alignments"]:
                        d[data[0]]["alignments"][data[5]] = [
                            (int(data[2]), int(data[3]))
                        ]
                    else:
                        d[data[0]]["alignments"][data[5]].append(
                            (int(data[2]), int(data[3]))
                        )
        # now we can go through each data set and calc coverage
        for contig, results in d.items():
            r = {}
            for hit, coords in results["alignments"].items():
                align_cov = 0
                for i in merge_coordinates(coords):
                    align_cov += i[1] - i[0]
                pct_cov = align_cov / results["length"]
                if pct_cov >= min_cov:
                    r[hit] = pct_cov
            if r:
                f[contig] = r
    else:
        f = None
    return f, d
