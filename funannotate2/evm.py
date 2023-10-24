#!/usr/bin/env python

import sys
import multiprocessing
import subprocess
import os
import time
import shutil
import argparse
import math
import pyfastx
from natsort import natsorted
from .interlap import InterLap
from .utilities import readBlocks, execute
from collections import defaultdict
from interlap import InterLap

from gfftk.gff import gff2dict, dict2gff3, dict2gff3alignments
from gfftk.fasta import softwrap
from collections import ChainMap
import bisect
from .utilities import runSubprocess, runProcessJob


def gene_blocks_to_interlap(input):
    # try to read EVM data as blocks and store in interlap by single interval
    inter = defaultdict(InterLap)
    with open(input, "r") as infile:
        for gene_model in readBlocks(infile, "\n"):
            # will be list of lines, so need to find gene line
            if len(gene_model) > 1:  # last line is has only a single element
                if gene_model[0] == "\n":
                    cols = gene_model[1].split("\t")
                else:
                    cols = gene_model[0].split("\t")
                inter[cols[0]].add((int(cols[3]), int(cols[4]), gene_model))
    return inter


def exonerate_blocks_to_interlap(input):
    # try to read EVM data as blocks and store in interlap by single interval
    inter = defaultdict(InterLap)
    with open(input, "r") as infile:
        for gene_model in readBlocks(infile, "\n"):
            coords = []
            if len(gene_model) < 2:
                continue
            for x in gene_model:
                if x == "\n":
                    continue
                else:
                    cols = x.split("\t")
                    coords.append(int(cols[3]))
                    coords.append(int(cols[4]))
            inter[cols[0]].add((min(coords), max(coords), gene_model))
    return inter


def blocks_to_interlap(input):
    # transcript alignments are not currently \n separated....
    # so read once to group and then build InterLap
    Results = {}
    with open(input, "r") as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            ID = None
            if ";" in cols[8]:
                info = cols[8].split(";")
            else:
                info = [cols[8]]
            for x in info:
                if x.startswith("ID="):
                    ID = x.replace("ID=", "")
            if not ID:
                continue
            if not ID in Results:
                Results[ID] = {
                    "coords": [int(cols[3]), int(cols[4])],
                    "raw": ["\n", line],
                    "contig": cols[0],
                }
            else:
                Results[ID]["raw"].append(line)
                Results[ID]["coords"].append(int(cols[3]))
                Results[ID]["coords"].append(int(cols[4]))
    # now build interlap object
    inter = defaultdict(InterLap)
    for k, v in sorted(Results.items()):
        inter[v["contig"]].add((min(v["coords"]), max(v["coords"]), v["raw"]))
    return inter


def worker(inputList):
    output = inputList[-2]
    logfile = inputList[-1]
    cmd = inputList[:-2]
    with open(output, "w") as outfile:
        with open(logfile, "w") as log:
            subprocess.call(cmd, stdout=outfile, stderr=log)


def safe_run(*args, **kwargs):
    """Call run(), catch exceptions."""
    try:
        worker(*args, **kwargs)
    except Exception as e:
        print(("error: %s run(*%r, **%r)" % (e, args, kwargs)))


def solve_partitions(values, interval=2000):
    # try to find best solution to paritioning
    # values is a list of tuples (start, stop, diff from prev, num of genes)
    result = []
    for i, x in enumerate(values):
        print(x)
        if x[2] >= interval:
            if not result:
                start = 1
            else:
                start = values[i - 1][1] + 1
            end = x[1] - 1
            if not result:
                result.append((start, end, x[2], x[3]))
            else:
                result.append((start, end, x[2], x[3] - result[-1][-1]))
    print("--------------")
    for y in result:
        print(y)


def create_partitions(
    fasta,
    genes,
    partition_list,
    proteins=False,
    transcripts=False,
    repeats=False,
    num=50,
    tmpdir=".",
    interval=2000,
    partitions=True,
    debug=False,
):
    # function to create EVM partition intervals that do not split genes
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)
    f_idx = fasta + ".idx"
    # SeqRecords = SeqIO.index_db(f_idx, fasta, 'fasta')
    SeqRecords = SeqIO.index(fasta, "fasta")
    PID = os.getpid()
    bedGenes = os.path.join(tmpdir, "genes.{}.bed".format(PID))
    superGenes = os.path.join(tmpdir, "genes.{}.supergenes.bed".format(PID))
    interGenes = gene_blocks_to_interlap(genes)
    if proteins:
        interProteins = exonerate_blocks_to_interlap(proteins)
    if transcripts:
        interTranscripts = blocks_to_interlap(transcripts)
    if repeats:
        interRepeats = blocks_to_interlap(repeats)
    Results = []
    with open(genes, "r") as infile:
        for line in infile:
            if line.startswith("#") or line.startswith("\n"):
                continue
            line = line.rstrip()
            cols = line.split("\t")
            if cols[2] == "gene":
                Results.append(
                    [cols[0], int(cols[3]), int(cols[4]), cols[8], cols[5], cols[6]]
                )
    # sort the results by contig and position
    ChrGeneCounts = {}
    totalGeneCount = 0
    sortedResults = natsorted(Results, key=lambda x: (x[0], x[1]))
    with open(bedGenes, "w") as outfile:
        for x in sortedResults:
            totalGeneCount += 1
            outfile.write(
                "{}\t{}\t{}\t{}\t{}\t{}\n".format(x[0], x[1], x[2], x[3], x[4], x[5])
            )
            if not x[0] in ChrGeneCounts:
                ChrGeneCounts[x[0]] = 1
            else:
                ChrGeneCounts[x[0]] += 1
    ChrNoGenes = len(SeqRecords) - len(ChrGeneCounts)
    superGeneCount = 0
    countbyContig = {}
    log.debug(
        "{:,} total contigs; skipping {:,} contigs with no genes".format(
            len(SeqRecords), ChrNoGenes
        )
    )
    if partitions:
        # now merge overlaping genes [strand] to get conservative locus boundaries
        cmd = ["bedtools", "merge", "-s", "-i", bedGenes]
        merged = {}
        with open(superGenes, "w") as outfile:
            for line in execute(cmd):
                superGeneCount += 1
                line = line.rstrip()
                if line.count("\t") != 2:
                    log.debug("Error parsing bedtools merge line:\n{}".format(line))
                    continue
                chr, start, end = line.split("\t")
                if not chr in countbyContig:
                    countbyContig[chr] = 1
                else:
                    countbyContig[chr] += 1
                outfile.write(
                    "{}\t{}\t{}\tSuperGene_{}\n".format(chr, start, end, superGeneCount)
                )
                if chr not in merged:
                    merged[chr] = [(int(start), int(end), -1, 1)]
                else:
                    diff = int(start) - merged[chr][-1][1]
                    merged[chr].append((int(start), int(end), diff, countbyContig[chr]))
        log.debug(
            "Merged {} genes into {} supergenes with bedtools".format(
                totalGeneCount, superGeneCount
            )
        )
        # parse Results and get coordinates to partitions
        Partitions = {}
        Commands = {}
        for k, v in natsorted(merged.items()):
            if not k in ChrGeneCounts:  # no genes, so can safely skip
                continue
            Partitions[k] = []
            next_start = None
            if len(v) > num:
                if not any(
                    z >= interval for z in [w[2] for w in v]
                ):  # this means there no intervals on this contig
                    Commands[k] = {"n": len(v)}
                    log.debug("{} --> {} bp".format(k, len(SeqRecords[k])))
                else:
                    chunks = math.ceil(len(v) / num)
                    num_genes = int(round(len(v) / chunks))
                    chunks = int(chunks)
                    for i in range(chunks):
                        if k in Commands:
                            continue
                        i = i + 1
                        if i == 1:
                            start = 1
                        else:
                            start = next_start
                        loc = i * num_genes
                        if i == chunks:
                            end = len(SeqRecords[k])
                        else:
                            if loc >= len(v):
                                end = len(SeqRecords[k])
                            else:
                                # new function to return the tuple to split
                                splitTup, idx = getBreakPoint(
                                    v, loc, direction="reverse", gap=interval
                                )
                                if not splitTup:
                                    splitTup, idx = getBreakPoint(
                                        v, loc, direction="forward", gap=interval
                                    )
                                end = splitTup[0] - 1
                                next_start = v[idx - 1][1] + 1
                        if not end:
                            Commands[k] = {"n": len(v)}
                            log.debug("{} --> {} bp".format(k, len(SeqRecords[k])))
                        else:
                            partLen = end - start
                            if partLen < 10000:
                                continue
                            Partitions[k].append((start, end))
                            partName = "{}_{}-{}".format(k, start, end)
                            Commands[partName] = {"n": num_genes, "chr": k}
                            log.debug("{} --> {} bp".format(partName, partLen))
                        start, end = (None,) * 2
            else:
                Commands[k] = {"n": len(v)}
                log.debug("{} --> {} bp".format(k, len(SeqRecords[k])))
        # now loop through partitions and write files for EVM
        with open(partition_list, "w") as partout:
            for chr, p in natsorted(Partitions.items()):
                chrDir = os.path.join(tmpdir, chr)
                if not os.path.isdir(chrDir):
                    os.makedirs(chrDir)
                if len(p) == 0:
                    partout.write(
                        "{}\t{}\t{}\n".format(chr, os.path.abspath(chrDir), "N")
                    )
                    chrFasta = os.path.join(chrDir, os.path.basename(fasta))
                    with open(chrFasta, "w") as fastaout:
                        fastaout.write(
                            ">{}\n{}\n".format(chr, softwrap(str(SeqRecords[chr].seq)))
                        )
                    genePred = os.path.join(chrDir, os.path.basename(genes))
                    RangeFinder(interGenes, chr, 1, len(SeqRecords[chr]), genePred)
                    if proteins:
                        protPred = os.path.join(chrDir, os.path.basename(proteins))
                        RangeFinder(
                            interProteins, chr, 1, len(SeqRecords[chr]), protPred
                        )
                    if transcripts:
                        tranPred = os.path.join(chrDir, os.path.basename(transcripts))
                        RangeFinder(
                            interTranscripts, chr, 1, len(SeqRecords[chr]), tranPred
                        )
                    if repeats:
                        repPred = os.path.join(chrDir, os.path.basename(repeats))
                        RangeFinder(interRepeats, chr, 1, len(SeqRecords[chr]), repPred)
                else:
                    for coords in p:
                        partDir = os.path.join(
                            chrDir, "{}_{}-{}".format(chr, coords[0], coords[1])
                        )
                        if not os.path.isdir(partDir):
                            os.makedirs(partDir)
                        partout.write(
                            "{}\t{}\t{}\t{}\n".format(
                                chr,
                                os.path.abspath(chrDir),
                                "Y",
                                os.path.abspath(partDir),
                            )
                        )
                        partFasta = os.path.join(partDir, os.path.basename(fasta))
                        with open(partFasta, "w") as fastaout:
                            fastaout.write(
                                ">{}\n{}\n".format(
                                    chr,
                                    softwrap(
                                        str(
                                            SeqRecords[chr].seq[
                                                coords[0] - 1 : coords[1]
                                            ]
                                        )
                                    ),
                                )
                            )
                        # split genes GFF3
                        genePred = os.path.join(partDir, "gene_predictions.gff3")
                        RangeFinder(interGenes, chr, coords[0], coords[1], genePred)
                        if proteins:
                            protPred = os.path.join(partDir, os.path.basename(proteins))
                            RangeFinder(
                                interProteins, chr, coords[0], coords[1], protPred
                            )
                        if transcripts:
                            tranPred = os.path.join(
                                partDir, os.path.basename(transcripts)
                            )
                            RangeFinder(
                                interTranscripts, chr, coords[0], coords[1], tranPred
                            )
                        if repeats:
                            repPred = os.path.join(partDir, os.path.basename(repeats))
                            RangeFinder(
                                interRepeats, chr, coords[0], coords[1], repPred
                            )
    else:
        Commands = {}
        with open(partition_list, "w") as partout:
            for chr in SeqRecords:
                if not chr in ChrGeneCounts:  # no genes so skip
                    continue
                Commands[chr] = {"n": len(SeqRecords[chr])}
                chrDir = os.path.join(tmpdir, chr)
                if not os.path.isdir(chrDir):
                    os.makedirs(chrDir)
                partout.write("{}\t{}\t{}\n".format(chr, os.path.abspath(chrDir), "N"))
                chrFasta = os.path.join(chrDir, os.path.basename(fasta))
                with open(chrFasta, "w") as fastaout:
                    fastaout.write(
                        ">{}\n{}\n".format(chr, softwrap(str(SeqRecords[chr].seq)))
                    )
                genePred = os.path.join(chrDir, os.path.basename(genes))
                RangeFinder(interGenes, chr, 1, len(SeqRecords[chr]), genePred)
                if proteins:
                    protPred = os.path.join(chrDir, os.path.basename(proteins))
                    RangeFinder(interProteins, chr, 1, len(SeqRecords[chr]), protPred)
                if transcripts:
                    tranPred = os.path.join(chrDir, os.path.basename(transcripts))
                    RangeFinder(
                        interTranscripts, chr, 1, len(SeqRecords[chr]), tranPred
                    )
                if repeats:
                    repPred = os.path.join(chrDir, os.path.basename(repeats))
                    RangeFinder(interRepeats, chr, 1, len(SeqRecords[chr]), repPred)
    SeqRecords.close()
    return Commands


def RangeFinder(input, chr, start, end, output, EVM=False):
    """
    if not EVM:
        EVM = os.environ['EVM_HOME']
    RFScript = os.path.join(EVM, 'EvmUtils', 'gff_range_retriever.pl')
    cmd = [RFScript, chr, str(start), str(end), 'ADJUST_TO_ONE']
    with open(output, 'w') as outfile:
        with open(os.path.abspath(input)) as infile:
            p = subprocess.Popen(cmd, stdin=infile, stdout=outfile)
            p.wait()
            outfile.flush()
    """
    hits = list(input[chr].find((start, end)))
    sortedHits = sorted(hits, key=lambda x: x[0])
    with open(output, "w") as outfile:
        for x in sortedHits:
            if x[0] < start or x[1] > end:
                continue
            adjust_coord = start - 1
            for line in x[2]:
                if line.startswith("#") or line.startswith("\n"):
                    outfile.write(line)
                else:
                    cols = line.split("\t")
                    cols[3] = "{}".format(int(cols[3]) - adjust_coord)
                    cols[4] = "{}".format(int(cols[4]) - adjust_coord)
                    outfile.write("{}".format("\t".join(cols)))


def getBreakPoint(tupList, idx, direction="reverse", gap=2000):
    # takes list of tuples of coords and a starting index (idx). finds closest
    # break point in between tuple coordSorted
    solution = False
    while not solution:
        try:
            start, end, diff, num_genes = tupList[idx]
        except IndexError:
            return False, idx
        if diff >= gap:
            # phase = int(round(diff/2))
            solution = tupList[idx]
            # print(idx, tupList[idx])
        else:
            if direction == "reverse":
                idx -= 1
            else:
                idx += 1
    return solution, idx


def main():
    # initialize script, log system info and cmd issue at runtime
    PID = os.getpid()
    if args.logfile:
        log_name = args.logfile
    else:
        log_name = "funannotate-evm.{}.log".format(PID)
    if os.path.isfile(log_name):
        os.remove(log_name)
    setupLogging(log_name)
    FNULL = open(os.devnull, "w")
    cmd_args = " ".join(sys.argv) + "\n"
    log.debug(cmd_args)

    # create output directory
    tmpdir = args.dir
    if tmpdir != ".":
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)
        os.makedirs(tmpdir)

    # set some EVM script locations
    perl = "perl"
    if args.EVM_HOME:
        EVM = args.EVM_HOME
    else:
        try:
            EVM = os.environ["EVM_HOME"]
        except:
            log.error("Could not find EVM_HOME environmental variable")
            sys.exit(1)

    Combine = os.path.join(EVM, "EvmUtils", "recombine_EVM_partial_outputs.pl")
    Convert = os.path.join(EVM, "EvmUtils", "convert_EVM_outputs_to_GFF3.pl")

    # split partitions
    partitions = os.path.join(tmpdir, "partitions_list.out")
    if args.no_partitions:
        log.info(
            "EVM: partitioning input to ~ {} genes per partition using min {} bp interval".format(
                args.gene_partition, args.interval
            )
        )
    else:
        log.info("EVM: partitioning each contig separately")
    cmdinfo = create_partitions(
        args.fasta,
        args.genes,
        partitions,
        proteins=args.proteins,
        transcripts=args.transcripts,
        repeats=args.repeats,
        tmpdir=tmpdir,
        num=args.gene_partition,
        interval=args.interval,
        partitions=args.no_partitions,
    )
    log.debug("Finished partitioning, generating command list")
    # sort the cmdinfo by number of putative genes and distribute into sub files
    num_workers = args.cpus - 1
    if num_workers < 1:
        num_workers = 1
    c = 0
    file_list = []
    tasks = 0
    for s in sorted(cmdinfo.items(), key=lambda x: x[1]["n"], reverse=True):
        key, d = s
        if "chr" in d:
            outputDir = os.path.abspath(os.path.join(tmpdir, d["chr"], key))
        else:
            outputDir = os.path.abspath(os.path.join(tmpdir, key))
        cmd = [
            os.path.join(EVM, "evidence_modeler.pl"),
            "-G",
            os.path.join(outputDir, os.path.basename(args.fasta)),
            "-g",
            os.path.join(outputDir, os.path.basename(args.genes)),
            "-w",
            os.path.abspath(args.weights),
            "--min_intron_length",
            str(args.min_intron),
            "--exec_dir",
            outputDir,
        ]
        if args.proteins:
            cmd += ["-p", os.path.join(outputDir, os.path.basename(args.proteins))]
        if args.transcripts:
            cmd += ["-e", os.path.join(outputDir, os.path.basename(args.transcripts))]
        if args.repeats:
            cmd += [
                "--repeats",
                os.path.join(outputDir, os.path.basename(args.repeats)),
            ]
        cmd += [
            os.path.join(outputDir, "evm.out"),
            os.path.join(outputDir, "evm.out.log"),
        ]
        file_list.append(cmd)

    # run runMultiProgress
    runMultiProgress(safe_run, file_list, num_workers, progress=args.progress)

    # now combine the paritions
    cmd4 = [
        perl,
        Combine,
        "--partitions",
        os.path.basename(partitions),
        "--output_file_name",
        "evm.out",
    ]
    runSubprocess(cmd4, tmpdir, log)

    # now convert to GFF3
    cmd5 = [
        perl,
        Convert,
        "--partitions",
        os.path.basename(partitions),
        "--output",
        "evm.out",
        "--genome",
        os.path.abspath(args.fasta),
    ]
    runSubprocess(cmd5, tmpdir, log)

    # now concatenate all GFF3 files together for a genome then
    log.info("Converting to GFF3 and collecting all EVM results")
    with open(args.out, "w") as out:
        for root, dirs, files in os.walk(tmpdir):
            for file in files:
                if file == "evm.out.gff3":
                    filename = os.path.join(root, file)
                    with open(filename, "r") as readfile:
                        shutil.copyfileobj(readfile, out)


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
        if not info["contig"] in all_genes:
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
    # run evidence modeler v2.1.1
    # check for EVM_HOME
    try:
        EVM = os.environ["EVM_HOME"]
    except:
        log.error("Could not find EVM_HOME environmental variable")
        raise SystemExit(1)

    # EvmUtils scripts
    Combine = os.path.join(EVM, "EvmUtils", "recombine_EVM_partial_outputs.pl")
    Convert = os.path.join(EVM, "EvmUtils", "convert_EVM_outputs_to_GFF3.pl")
    RunEVM = os.path.join(EVM, "EvmUtils", "evidence_modeler.pl")

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
                    end_position = int(
                        (merged_loci[end_idx + 1][0] + merged_loci[end_idx][1]) / 2
                    )
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
    subprocess.call(
        [
            Combine,
            "--partitions",
            os.path.join(tmpdir, "partitions_list.out"),
            "--output_file_name",
            "evm.out",
        ]
    )

    # convert to gff3
    subprocess.call(
        [
            Convert,
            "--partitions",
            os.path.join(tmpdir, "partitions_list.out"),
            "--output_file_name",
            "evm.out",
            "--genome",
            genome,
        ]
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
                partout.write(
                    "{}\t{}\t{}\n".format(contig, os.path.abspath(chrDir), "N")
                )
                # write the fasta file
                chrFasta = os.path.join(chrDir, "genome.fasta")
                with open(chrFasta, "w") as fastaout:
                    fastaout.write(
                        ">{}\n{}\n".format(contig, softwrap(Genome[contig].seq))
                    )
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
                        output=tranPred,
                        alignments="protein",
                        newline=True,
                    )
                    subcmd += ["--protein_alignments", os.path.basename(protPred)]

                if len(v["repeats"]) > 0:
                    print("more to do here")

            elif n_parts[contig] > 1:  # chr is split into chunks
                partDir = os.path.join(chrDir, k)
                if not os.path.isdir(partDir):
                    os.makedirs(partDir)
                partout.write(
                    "{}\t{}\t{}\t{}\n".format(
                        contig, os.path.abspath(chrDir), "Y", os.path.abspath(partDir)
                    )
                )
                part_interval = (
                    int(k.split(":")[1].split("-")[0]),
                    int(k.split(":")[1].split("-")[1]),
                )
                # write the fasta file
                chrFasta = os.path.join(partDir, "genome.fasta")
                with open(chrFasta, "w") as fastaout:
                    fastaout.write(
                        ">{}\n{}\n".format(
                            contig, softwrap(Genome.fetch(contig, part_interval))
                        )
                    )
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
                        v["protein_evidence"], output=tranPred, alignments="protein"
                    )
                    subcmd += ["--protein_alignments", os.path.basename(protPred)]

                if len(v["repeats"]) > 0:
                    print("more to do here")
            cmds.append(subcmd)
    return cmds
