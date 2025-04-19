#!/usr/bin/env python

import bisect
import math
import os
import shutil
import subprocess
import sys
from collections import ChainMap, defaultdict

import pyfastx
from gfftk.fasta import softwrap
from gfftk.gff import dict2gff3, dict2gff3alignments, gff2dict

from .interlap import InterLap
from .utilities import checkfile, runProcessJob, runSubprocess


def merge(coords):
    saved = list(coords[0])
    for st, en in sorted([sorted(t) for t in coords]):
        if st <= saved[1]:
            saved[1] = max(saved[1], en)
        else:
            yield tuple(saved)
            saved[0] = st
            saved[1] = en
    yield tuple(saved)


def get_loci(annot_dict):
    # input is annotation dictionary
    # so loop through and add to contig keyed dictionary
    all_genes = {}
    for gene, info in annot_dict.items():
        if info["contig"] not in all_genes:
            all_genes[info["contig"]] = [info["location"]]
        else:
            all_genes[info["contig"]].append(info["location"])
    # now we need to merge overlapping to generate loci
    loci = {}
    for chr, g in all_genes.items():
        merged = list(merge(g))
        loci[chr] = merged
    return loci


def run_evm_partition(cmd_list, **kwargs):
    idx = cmd_list.index("--exec_dir") + 1
    outevm = os.path.join(cmd_list[idx], "evm.out")
    logevm = os.path.join(cmd_list[idx], "evm.log")
    with open(outevm, "w") as outfile:
        with open(logevm, "w") as logfile:
            proc = subprocess.Popen(cmd_list, stdout=outfile, stderr=logfile, **kwargs)
            stdout, stderr = proc.communicate()
            if proc.returncode != 0:
                print("CMD ERROR: {}".format(" ".join(cmd_list)))
                raise SystemExit(1)


def evm_consensus(
    genome,
    gene_preds,
    prot_aligns,
    transcript_aligns,
    weights,
    consensus,
    n_genes=35,
    intergenic=1500,
    tmpdir="/tmp",
    log=sys.stderr.write,
    cpus=1,
):
    # run evidence modeler v2.1.1 or older version, basically need to see where evidence_modeler.pl is located
    # check for EVM_HOME
    try:
        EVM = os.environ["EVM_HOME"]
    except KeyError:
        log.error("Could not find EVM_HOME environmental variable")
        raise SystemExit(1)

    # EvmUtils scripts
    Combine = os.path.join(EVM, "EvmUtils", "recombine_EVM_partial_outputs.pl")
    Convert = os.path.join(EVM, "EvmUtils", "convert_EVM_outputs_to_GFF3.pl")
    RunEVM = os.path.join(EVM, "EvmUtils", "evidence_modeler.pl")
    if not checkfile(RunEVM):
        RunEVM = os.path.join(EVM, "evidence_modeler.pl")
        if not checkfile(RunEVM):
            log.error(
                f"EVM installation problem, unable to location evidence_modeler.pl script in {EVM}"
            )

    # load the genome
    fa = pyfastx.Fasta(genome)

    # parse the gene_preds into loci so you know where to partition
    Preds = []
    n_models = 0
    sources = {"ABINITIO_PREDICTION": set(), "PROTEIN": set(), "TRANSCRIPT": set()}
    log.info(gene_preds)
    for g in gene_preds:
        parsed_genes = gff2dict(os.path.abspath(g), genome)
        Preds.append(parsed_genes)
        n_models += len(parsed_genes)
        for k, v in parsed_genes.items():
            sources["ABINITIO_PREDICTION"].add(v["source"])
    Proteins = {}
    log.info(prot_aligns)
    if prot_aligns:
        for p in prot_aligns:
            Proteins = gff2dict(os.path.abspath(p), genome, gff_format="alignment")
            for k, v in Proteins.items():
                sources["PROTEIN"].add(v["source"])
    Transcripts = {}
    log.info(transcript_aligns)
    if transcript_aligns:
        for t in transcript_aligns:
            Transcripts = gff2dict(os.path.abspath(t), genome, gff_format="alignment")
            for k, v in Transcripts.items():
                sources["TRANSCRIPT"].add(v["source"])
    # write weights file
    weightsfile = os.path.join(tmpdir, "weights.txt")
    w = {}
    for x in weights:
        tool, wt = x.split(":", 1)
        w[tool] = wt
    with open(weightsfile, "w") as wfile:
        for k, v in sources.items():
            for t in v:
                wt_score = w.get(t, 1)
                wfile.write("{}\t{}\t{}\n".format(k, t, wt_score))

    # merge parsed predictions
    log.info(f"Merging gene predictions from {len(Preds)} source files\n{sources}")
    Genes = dict(ChainMap(*Preds))
    # now build loci from gene models, returns loci[chr] = {"+": plus_bins, "-": minus_bins}
    loci = get_loci(Genes)
    Partitions = {}
    for k, merged_loci in loci.items():
        if len(merged_loci) < n_genes * 2:
            Partitions[k] = [(1, len(fa[k]))]
        else:
            splits = []
            chunks = int(math.ceil(len(merged_loci) / n_genes))
            steps = int(len(merged_loci) / chunks)
            step_values = [steps * x for x in range(1, chunks)]
            possible_split_pos = []
            for i in range(1, len(merged_loci)):
                intgen = merged_loci[i][0] - merged_loci[i - 1][1]
                if intgen >= intergenic:
                    possible_split_pos.append(i)
            for i, z in enumerate(step_values):
                idx = bisect.bisect_left(possible_split_pos, z)
                end_idx = possible_split_pos[idx]
                if merged_loci[-1] == merged_loci[end_idx]:  # then this is the last one
                    start = splits[-1][1] + 1
                    end_position = len(fa[k])
                else:
                    end_position = int((merged_loci[end_idx + 1][0] + merged_loci[end_idx][1]) / 2)
                    if i == 0:
                        start = 1
                    else:
                        start = splits[-1][1] + 1
                        if i == len(step_values) - 1:
                            end_position = len(fa[k])
                splits.append((start, end_position))
            Partitions[k] = splits
    # now split data and generate commands
    cmdlist = cmd_splitter(
        Partitions, fa, Genes, Transcripts, Proteins, {}, tmpdir, weightsfile, RunEVM
    )
    # run these in multi process
    # runProcessJob(function, inputList, cpus=2)
    log.info(f"EVM jobs running in {tmpdir}")
    runProcessJob(run_evm_partition, cmdlist, cpus=cpus)

    # now combine the outputs from partitions
    runSubprocess(
        [
            Combine,
            "--partitions",
            os.path.join(tmpdir, "partitions_list.out"),
            "--output_file_name",
            "evm.out",
        ],
        log,
        only_failed=True,
    )
    # convert to gff3
    runSubprocess(
        [
            Convert,
            "--partitions",
            os.path.join(tmpdir, "partitions_list.out"),
            "--output_file_name",
            "evm.out",
            "--genome",
            genome,
        ],
        log,
        only_failed=True,
    )

    # gather up all gff3 results
    with open(consensus, "w") as out:
        for root, dirs, files in os.walk(tmpdir):
            for file in files:
                if file == "evm.out.gff3":
                    filename = os.path.join(root, file)
                    with open(filename, "r") as readfile:
                        shutil.copyfileobj(readfile, out)


def adjust_coords(m, offset):
    # take a gene model dictionary and adjust coordinates
    new = m.copy()
    new["location"] = ((m["location"][0] - offset), (m["location"][1] - offset))
    if "mRNA" in new:
        n = []
        for x in m["mRNA"]:
            n.append([(t[0] - offset, t[1] - offset) for t in x])
        new["mRNA"] = n
    if "CDS" in new:
        n = []
        for x in m["CDS"]:
            n.append([(t[0] - offset, t[1] - offset) for t in x])
        new["CDS"] = n
    if "5UTR" in new:
        n = []
        for x in m["5UTR"]:
            n.append([(t[0] - offset, t[1] - offset) for t in x])
        new["5UTR"] = n
    if "3UTR" in new:
        n = []
        for x in m["3UTR"]:
            n.append([(t[0] - offset, t[1] - offset) for t in x])
        new["3UTR"] = n
    return new


def cmd_splitter(
    Partitions,
    Genome,
    Genes,
    Transcripts,
    Proteins,
    Repeats,
    writedir,
    weights_file,
    evmcmd,
    min_intron=10,
):
    # data in dicts, split by location and write files to their loc and write evm commands
    # lets split the data for each partition and then write
    cmds = []
    data = {}
    p_inter = defaultdict(InterLap)
    n_parts = {}
    for k, v in Partitions.items():
        n_parts[k] = len(v)
        for c in v:
            part_name = f"{k}:{c[0]}-{c[1]}"
            p_inter[k].add((c[0], c[1], part_name))
            data[part_name] = {
                "genes": {},
                "protein_evidence": {},
                "transcript_evidence": {},
                "repeats": {},
            }
    # now we can go through the models and add to each dictionary
    for k, v in Genes.items():
        hits = list(p_inter[v["contig"]].find(v["location"]))
        if len(hits) == 1:
            part = hits[0][2]
            offset = int(part.split(":")[1].split("-")[0]) - 1
            if offset == 0:  # then do nothing
                data[part]["genes"][k] = v
            elif offset > 0:  # then we have to adjust coordinates
                data[part]["genes"][k] = adjust_coords(v, offset)
        else:
            print(
                "This should not happen, {} {} is not in the interlap object".format(
                    v["contig"], v["location"]
                )
            )
            raise SystemExit(1)
    # now we can go through the evidence
    for k, v in Transcripts.items():
        hits = list(p_inter[v["contig"]].find(v["location"]))
        if len(hits) > 0:
            for h in hits:
                part = h[2]
                offset = int(part.split(":")[1].split("-")[0]) - 1
                if offset == 0:  # then do nothing
                    data[part]["transcript_evidence"][k] = v
                elif offset > 0:  # then we have to adjust coordinates
                    data[part]["transcript_evidence"][k] = adjust_coords(v, offset)
        else:
            print(
                "This should not happen, {} {} is not in the interlap object".format(
                    v["contig"], v["location"]
                )
            )
            raise SystemExit(1)
    for k, v in Proteins.items():
        hits = list(p_inter[v["contig"]].find(v["location"]))
        if len(hits) > 0:
            for h in hits:
                part = h[2]
                offset = int(part.split(":")[1].split("-")[0]) - 1
                if offset == 0:  # then do nothing
                    data[part]["protein_evidence"][k] = v
                elif offset > 0:  # then we have to adjust coordinates
                    data[part]["protein_evidence"][k] = adjust_coords(v, offset)
        else:
            print(
                "This should not happen, {} {} is not in the interlap object".format(
                    v["contig"], v["location"]
                )
            )
            raise SystemExit(1)
    for k, v in Repeats.items():
        hits = list(p_inter[v["contig"]].find(v["location"]))
        if len(hits) > 0:
            for h in hits:
                part = h[2]
                offset = int(part.split(":")[1].split("-")[0]) - 1
                if offset == 0:  # then do nothing
                    data[part]["repeats"][k] = v
                elif offset > 0:  # then we have to adjust coordinates
                    data[part]["repeats"][k] = adjust_coords(v, offset)
                data[part]["repeats"][k] = v
        else:
            print(
                "This should not happen, {} {} is not in the interlap object".format(
                    v["contig"], v["location"]
                )
            )
            raise SystemExit(1)
    # write the data and generate the evidence modeler commands
    # this is tab delimited with following values:
    # ($accession, $base_dir, $partitioned, $partition_dir)
    partitions_file = os.path.join(writedir, "partitions_list.out")
    with open(partitions_file, "w") as partout:
        for k, v in data.items():
            contig = k.split(":")[0]
            chrDir = os.path.join(writedir, contig)
            if not os.path.isdir(chrDir):
                os.makedirs(chrDir)
            # write entry in partitons list and the components
            subcmd = [
                evmcmd,
                "--weights",
                weights_file,
                "--min_intron_length",
                str(min_intron),
            ]
            if n_parts[contig] == 1:  # entire chromosome as parition
                partout.write("{}\t{}\t{}\n".format(contig, os.path.abspath(chrDir), "N"))
                # write the fasta file
                chrFasta = os.path.join(chrDir, "genome.fasta")
                with open(chrFasta, "w") as fastaout:
                    fastaout.write(">{}\n{}\n".format(contig, softwrap(Genome[contig].seq)))
                # write the gene predictions
                genePred = os.path.join(chrDir, "gene_predictions.gff3")
                dict2gff3(v["genes"], output=genePred, newline=True)

                # setup command, others could be extra
                subcmd += [
                    "--exec_dir",
                    chrDir,
                    "--genome",
                    os.path.basename(chrFasta),
                    "--gene_predictions",
                    os.path.basename(genePred),
                ]
                if len(v["transcript_evidence"]) > 0:
                    tranPred = os.path.join(chrDir, "transcript_alignments.gff3")
                    dict2gff3alignments(
                        v["transcript_evidence"],
                        output=tranPred,
                        alignments="transcript",
                        newline=True,
                    )
                    subcmd += ["--transcript_alignments", os.path.basename(tranPred)]
                if len(v["protein_evidence"]) > 0:
                    protPred = os.path.join(chrDir, "protein_alignments.gff3")
                    dict2gff3alignments(
                        v["protein_evidence"],
                        output=protPred,
                        alignments="protein",
                        newline=True,
                    )
                    subcmd += ["--protein_alignments", os.path.basename(protPred)]

                if len(v["repeats"]) > 0:
                    print("more to do here")

            elif n_parts[contig] > 1:  # chr is split into chunks
                part_interval = (
                    int(k.split(":")[1].split("-")[0]),
                    int(k.split(":")[1].split("-")[1]),
                )
                if part_interval[0] > part_interval[1]:
                    continue
                partDir = os.path.join(chrDir, k)
                if not os.path.isdir(partDir):
                    os.makedirs(partDir)
                partout.write(
                    "{}\t{}\t{}\t{}\n".format(
                        contig, os.path.abspath(chrDir), "Y", os.path.abspath(partDir)
                    )
                )

                # write the fasta file
                chrFasta = os.path.join(partDir, "genome.fasta")
                try:
                    with open(chrFasta, "w") as fastaout:
                        fastaout.write(
                            ">{}\n{}\n".format(
                                contig, softwrap(Genome.fetch(contig, part_interval))
                            )
                        )
                except ValueError:
                    print(f"ERROR fetching contig: {contig} {part_interval}")
                    raise SystemExit(1)
                # write the gene predictions
                genePred = os.path.join(partDir, "gene_predictions.gff3")
                dict2gff3(v["genes"], output=genePred)

                # setup command, others could be extra
                subcmd += [
                    "--exec_dir",
                    partDir,
                    "--genome",
                    os.path.basename(chrFasta),
                    "--gene_predictions",
                    os.path.basename(genePred),
                ]
                if len(v["transcript_evidence"]) > 0:
                    tranPred = os.path.join(partDir, "transcript_alignments.gff3")
                    dict2gff3alignments(
                        v["transcript_evidence"],
                        output=tranPred,
                        alignments="transcript",
                    )
                    subcmd += ["--transcript_alignments", os.path.basename(tranPred)]
                if len(v["protein_evidence"]) > 0:
                    protPred = os.path.join(partDir, "protein_alignments.gff3")
                    dict2gff3alignments(
                        v["protein_evidence"], output=protPred, alignments="protein"
                    )
                    subcmd += ["--protein_alignments", os.path.basename(protPred)]

                if len(v["repeats"]) > 0:
                    print("more to do here")
            cmds.append(subcmd)
    return cmds
