import sys
import subprocess
import uuid
import os
from gfftk.gff import gff2dict, dict2gff3, dict2gff3alignments
from gapmm2.align import splice_aligner, paf2gff3
from .utilities import runSubprocess


def split_evidence_and_genes(gff, fasta, evidence, genes, gff_format="alignment"):
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
    # function to run spliced transcript alignment using gapmm2 (mm2 + edlib refinement)
    # generate an EVM compatible GFF alignment file, parse data and pull out any full length gene models
    log.info("Aligning transcript evidence to the genome assembly with gapmm2")
    t_aligned, t_stats = splice_aligner(
        fasta, transcripts, threads=cpus, max_intron=max_intron
    )
    paf2gff3(t_aligned, evidence)
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
    runSubprocess(cmd, log, stdout=mini_tmp)
    # miniprot uses non-standard gff, so load with gfftk parse
    # now we need to parse this evidence GFF3 file and look for full length
    n_aligns, n_genes = split_evidence_and_genes(
        mini_tmp, fasta, evidence, genes, gff_format="miniprot"
    )
    log.info(f"Generated {n_aligns} alignments: {n_genes} were valid gene models")
    os.remove(mini_tmp)
