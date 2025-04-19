import errno
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import types
import uuid
from collections import OrderedDict

import numpy as np
import pyfastx
from gfftk.fasta import fasta2dict, softwrap
from gfftk.gff import dict2gff3, gff2dict, gtf2dict, validate_models
from natsort import natsorted

from .config import env
from .genbank import gff3_to_gbio
from .interlap import InterLap
from .utilities import (
    checkfile,
    execute,
    print_table,
    readBlocks,
    runSubprocess,
    which_path,
)


class reversor:
    def __init__(self, obj):
        self.obj = obj

    def __eq__(self, other):
        return other.obj == self.obj

    def __lt__(self, other):
        return other.obj < self.obj


def run_trnascan(genome, predictions, cpus=1, folder="/tmp", log=sys.stderr.write):
    """
    Run tRNAscan-SE to predict tRNA genes in a genome and convert the results to GFF3 format.

    This function processes the genome file to predict tRNA genes, filters predictions based on length, and writes the results in GFF3 format.

    Args:
        genome (str): Path to the input genome file.
        predictions (str): Path to save the predicted tRNA genes in GFF3 format.
        cpus (int, optional): Number of CPUs to use for prediction. Defaults to 1.
        folder (str, optional): Temporary folder to store intermediate files. Defaults to '/tmp'.
        log (function, optional): Logger function for debug and error messages. Defaults to sys.stderr.write.

    Returns:
        int: Returns 0 if successful.
    """

    def _trnaraw2gff3(raw, gff3):
        """
        Sequence   		tRNA   	Bounds 	tRNA	Anti	Intron Bounds	Inf
        Name       	tRNA #	Begin  	End    	Type	Codon	Begin	End	Score	Note
        --------   	------	-----  	------ 	----	-----	-----	----	------	------
        scaffold_1 	1	661826 	661896 	Gly	GCC	0	0	64.0
        scaffold_1 	2	1175981	1176054	Ile	AAT	0	0	74.0
        scaffold_1 	3	1306530	1306629	Ser	TGA	1306567	1306585	66.8
        scaffold_1 	4	1343311	1343240	iMet	CAT	0	0	69.3
        scaffold_1 	5	1331666	1331595	Gly	CCC	0	0	47.0
        scaffold_1 	6	1225093	1225001	Ser	CGA	1225056	1225045	73.4
        scaffold_1 	7	1078648	1078577	Ala	AGC	0	0	68.8
        """
        with open(gff3, "w") as outfile:
            outfile.write("##gff-version 3\n")
            with open(raw, "r") as infile:
                data = infile.readlines()
            tdata = [[x.strip() for x in d.rstrip().split("\t")] for d in data[3:]]
            for i, pred in enumerate(tdata):
                feature = "gene"
                (
                    contig,
                    n,
                    start,
                    end,
                    t_type,
                    anticodon,
                    intron_start,
                    intron_end,
                    score,
                ) = pred[:9]
                if len(pred) == 10:
                    if pred[-1] == "pseudo":
                        feature = "pseudogene"
                start = int(start)
                end = int(end)
                intron_start = int(intron_start)
                intron_end = int(intron_end)
                gene_name = f"{t_type}_{i + 1}"
                trna_name = f"{gene_name}_tRNA"
                product = f"tRNA-{t_type}"
                note = f"Predicted {anticodon} anticodon"
                if start < end:
                    strand = "+"
                    outfile.write(
                        f"{contig}\ttRNAScan-SE\t{feature}\t{start}\t{end}\t{score}\t{strand}\t.\tID={gene_name};\n"
                    )
                    outfile.write(
                        f"{contig}\ttRNAScan-SE\ttRNA\t{start}\t{end}\t{score}\t{strand}\t.\tID={trna_name};Parent={gene_name};product={product};note={note};\n"
                    )
                    if intron_start > 0 and intron_end > 0:
                        outfile.write(
                            f"{contig}\ttRNAScan-SE\texon\t{start}\t{intron_start - 1}\t.\t{strand}\t.\tID={trna_name}.exon1;Parent={trna_name};\n"
                        )
                        outfile.write(
                            f"{contig}\ttRNAScan-SE\texon\t{intron_end + 1}\t{end}\t.\t{strand}\t.\tID={trna_name}.exon2;Parent={trna_name};\n"
                        )
                    else:
                        outfile.write(
                            f"{contig}\ttRNAScan-SE\texon\t{start}\t{end}\t.\t{strand}\t.\tID={trna_name}.exon1;Parent={trna_name};\n"
                        )
                else:
                    strand = "-"
                    outfile.write(
                        f"{contig}\ttRNAScan-SE\t{feature}\t{end}\t{start}\t{score}\t{strand}\t.\tID={gene_name};\n"
                    )
                    outfile.write(
                        f"{contig}\ttRNAScan-SE\ttRNA\t{end}\t{start}\t{score}\t{strand}\t.\tID={trna_name};Parent={gene_name};product={product};note={note};\n"
                    )
                    if intron_start > 0 and intron_end > 0:
                        outfile.write(
                            f"{contig}\ttRNAScan-SE\texon\t{end}\t{intron_end + 1}\t.\t{strand}\t.\tID={trna_name}.exon1;Parent={trna_name};\n"
                        )
                        outfile.write(
                            f"{contig}\ttRNAScan-SE\texon\t{intron_start - 1}\t{start}\t.\t{strand}\t.\tID={trna_name}.exon2;Parent={trna_name};\n"
                        )
                    else:
                        outfile.write(
                            f"{contig}\ttRNAScan-SE\texon\t{end}\t{start}\t.\t{strand}\t.\tID={trna_name}.exon1;Parent={trna_name};\n"
                        )

    # generate training directory ontop of the dir that is passed
    tmpdir = os.path.join(folder, f"trnascan-{uuid.uuid4()}")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    os.makedirs(tmpdir)
    tRNAout = os.path.join(tmpdir, "tRNAscan.out")
    tRNAlenOut = os.path.join(tmpdir, "tRNAscan.len-filtered.out")
    cmd = [
        "tRNAscan-SE",
        "-o",
        os.path.abspath(tRNAout),
        "--thread",
        str(cpus),
        os.path.abspath(genome),
    ]
    with tempfile.TemporaryDirectory() as tmpdirname:
        runSubprocess(cmd, log, cwd=tmpdirname)

    # enforce NCBI length rules
    with open(tRNAlenOut, "w") as lenOut:
        with open(tRNAout, "r") as infile:
            for line in infile:
                if (
                    line.startswith("Sequence")
                    or line.startswith("Name")
                    or line.startswith("--------")
                ):
                    lenOut.write(f"{line}")
                else:
                    cols = line.split("\t")
                    start = cols[2]
                    end = cols[3]
                    if int(start) < int(end):
                        length = abs(int(end) - int(start))
                    else:
                        length = abs(int(start) - int(end))
                    if length < 50 or length > 150:
                        continue
                    else:
                        lenOut.write(f"{line}")
    # convert to GFF3 format
    _trnaraw2gff3(tRNAlenOut, predictions)
    # if os.path.isdir(tmpdir):
    #    shutil.rmtree(tmpdir)

    return 0


def train_glimmerhmm(
    genome,
    train_gff,
    test_models,
    species=None,
    minintron=10,
    maxintron=3000,
    folder="/tmp",
    asm_size=None,
    log=sys.stderr.write,
):
    """
    Train a GlimmerHMM model using the provided genome and training data, and optimize parameters through grid search.

    This function prepares the training data, runs the initial training, and performs parameter optimization to improve model accuracy.

    Args:
        genome (str): Path to the genome file.
        train_gff (str or dict): Training data in GFF format or a dictionary.
        test_models (dict): Dictionary of test models.
        species (str, optional): Species name. Defaults to None.
        minintron (int, optional): Minimum intron size. Defaults to 10.
        maxintron (int, optional): Maximum intron size. Defaults to 3000.
        folder (str, optional): Directory to store temporary files. Defaults to "/tmp".
        asm_size (float, optional): Assembly size for parameter calculation. Defaults to None.
        log (function, optional): Logger for messages. Defaults to sys.stderr.write.

    Returns:
        dict: A dictionary containing the training results and the location of the trained model.
    """

    def _dict2glimmer(input, exon_out, fasta_out):
        # take funannotate dictionary convert to glimmer training format
        with open(fasta_out, "w") as fastaout:
            with open(exon_out, "w") as outfile:
                for k, v in input.items():
                    try:
                        fastaout.write(f">{k}\n{softwrap(v['train_dna'])}\n")
                    except KeyError:
                        print(k, v)
                        raise SystemExit(1)
                    for c in v["train_coords"]:
                        if v["strand"] == "+":
                            outfile.write("{:} {:} {:}\n".format(k, c[0], c[1]))
                        else:
                            outfile.write("{:} {:} {:}\n".format(k, c[1], c[0]))
                    outfile.write("\n")

    def _glimmer_config(input):
        # parse glimmer config file and return a dictionary
        cfg = {}
        with open(input, "r") as infile:
            for line in infile:
                line = line.rstrip()
                if " " in line:
                    key, value = line.split(" ", 1)
                    if key not in cfg:
                        cfg[key] = value
        return cfg

    def _new_glimmer_train_data(init_train, cfg, new_train):
        # init_train is the existing directory
        # cfg is the glimmer config dictionary
        # new_train is the new directory to be created
        if not os.path.isdir(new_train):
            os.makedirs(new_train)
        with open(os.path.join(new_train, "train_0_100.cfg"), "w") as outfile:
            for k, v in cfg.items():
                outfile.write("{:} {:}\n".format(k, v))
        for f in os.listdir(init_train):
            if ".cfg" not in f:
                shutil.copy(os.path.join(init_train, f), new_train)

    def _load_training_parameters(train_dir):
        # load training parameters from glimmer config file
        cfg = _glimmer_config(os.path.join(train_dir, "train_0_100.cfg"))
        params = {
            "split_penalty": {0: {}, 1: {}, 2: {}, 3: {}, 4: {}},
            "intergenic_val": {},
            "intergenic_penalty": {0: {}, 10: {}, 25: {}, 50: {}, 100: {}},
            "MeanIntergen": {},
            "BoostExon": {0: {}, 5: {}, 10: {}, 12: {}, 15: {}},
            "BoostSplice": {0: {}, 5: {}, 10: {}, 12: {}, 15: {}},
            "BoostSgl": {0: {}, 5: {}, 10: {}, 12: {}, 15: {}},
            "use_intron_distrib": {0: {}, 1: {}},
        }
        # set the params based on the initial tranining in cfg
        iv = float(cfg["intergenic_val"])
        if iv <= 200:
            params["intergenic_val"] = {
                iv: {},
                iv + 200: {},
                iv + 500: {},
                iv + 1000: {},
                iv + 1500: {},
            }
        elif iv <= 500:
            params["intergenic_val"] = {
                iv - 200: {},
                iv: {},
                iv + 200: {},
                iv + 500: {},
                iv + 1000: {},
            }
        else:
            params["intergenic_val"] = {
                iv - 500: {},
                iv - 200: {},
                iv: {},
                iv + 200: {},
                iv + 500: {},
            }
        mi = float(cfg["MeanIntergen"])
        if mi <= 1000:
            params["MeanIntergen"] = {
                mi: {},
                mi + 500: {},
                mi + 1000: {},
                mi + 1500: {},
                mi + 2000: {},
            }
        elif mi <= 2000:
            params["MeanIntergen"] = {
                mi - 500: {},
                mi: {},
                mi + 500: {},
                mi + 1000: {},
                mi + 1500: {},
            }
        else:
            params["MeanIntergen"] = {
                mi - 1000: {},
                mi - 500: {},
                mi: {},
                mi + 500: {},
                mi + 1000: {},
            }
        return params

    # generate training directory ontop of the dir that is passed
    tmpdir = os.path.join(folder, "glimmerhmm")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    os.makedirs(tmpdir)

    # generate glimmer training input
    # load gff3 into dictionary
    if isinstance(train_gff, str):
        Genes = gff2dict(os.path.abspath(train_gff), os.path.abspath(genome), logger=log)
    elif isinstance(train_gff, dict):
        Genes = train_gff
    glimmExons = os.path.join(tmpdir, "glimmer.exons")
    glimmFasta = os.path.join(tmpdir, "glimmer.fasta")

    # calculate glimmer exon file needed for training
    _dict2glimmer(Genes, glimmExons, glimmFasta)

    TrainData = os.path.abspath(os.path.join(tmpdir, "train"))

    time_start = time.time()
    if not os.path.isdir(TrainData):
        # now run trainGlimmerHMM
        cmd = [
            "trainGlimmerHMM",
            os.path.abspath(glimmFasta),
            os.path.abspath(glimmExons),
            "-d",
            "train",
        ]
        runSubprocess(cmd, log, cwd=tmpdir, only_failed=True)
    train_end = time.time()
    train_elapsed = time.strftime("%H:%M:%S", time.gmtime(train_end - time_start))

    # now we can test the trained model
    test_training(TrainData, test_models, tmpdir, tool="glimmerhmm")
    # log.info(
    #    "Initial training completed in {}s\n{}".format(
    #        train_elapsed, json.dumps(init_train, indent=2, default=str)
    #    )
    # )

    # we can start with some obvious changes/updates to the config, add split penalty and ensure intergenic is larger than default
    glimmCFG = os.path.join(tmpdir, "train", "train_0_100.cfg")
    os.rename(glimmCFG, glimmCFG + ".bak")
    # calculate mean intergenic size based on genome, but never lower than 1 kb
    min_intergenic = 200.0
    mean_intergenic = 1000.0
    if asm_size:
        min_intergenic = round((asm_size / 1e6), 0)
        if min_intergenic < 200.0:
            min_intergenic = 200.0
        mean_intergenic = round((asm_size / 1e6) * 3, 0)
        if mean_intergenic < 1000.0:
            mean_intergenic = 1000.0
    with open(glimmCFG, "w") as outfile:
        with open(glimmCFG + ".bak", "r") as infile:
            for line in infile:
                if line.startswith("split_penalty"):
                    outfile.write("split_penalty 2\n")
                    outfile.write("use_intron_distrib 0\n")
                elif line.startswith("MeanIntergen"):
                    outfile.write(f"MeanIntergen {mean_intergenic}\n")
                    outfile.write(f"intergenic_val {min_intergenic}\n")
                    outfile.write("intergenic_penalty 25\n")
                elif line.startswith("BoostExon"):
                    outfile.write(line)
                    outfile.write("BoostSplice 0\n")
                    outfile.write("BoostSgl 0\n")
                else:
                    outfile.write(line)
    # write the test models to their own directory
    if not os.path.isdir(os.path.join(tmpdir, "test")):
        os.makedirs(os.path.join(tmpdir, "test"))
        for k, v in test_models.items():
            with open(os.path.join(tmpdir, "test", f"{k}.fa"), "w") as outfile:
                outfile.write(f">{k}\n{softwrap(v['train_dna'])}\n")

    # if optimize passed, then insert grid search logic here
    base_cfg = _glimmer_config(glimmCFG)
    optimize_dir = os.path.join(tmpdir, "optimize")
    if not os.path.isdir(optimize_dir):
        os.makedirs(optimize_dir)

    # get the parameters to try to optimize
    grid_search = _load_training_parameters(TrainData)
    # now run the grid search
    best_params = {}
    for k, v in grid_search.items():
        res = {}
        for i in v.keys():
            new_cfg = base_cfg.copy()
            new_cfg[k] = i
            new_train = os.path.join(optimize_dir, f"train_{k}_{i}")
            _new_glimmer_train_data(TrainData, new_cfg, new_train)
            r = test_training(new_train, test_models, tmpdir, tool="glimmerhmm")
            res[i] = r
        # print("grid search for {}".format(k))
        # print(json.dumps(res, indent=2))
        sRes = sorted(
            res.items(),
            key=lambda x: (
                x[1]["average_aed"],
                reversor(x[1]["nucleotide_precision"]),
                reversor(x[1]["exon_precision"]),
                reversor(x[1]["gene_precision"]),
            ),
        )
        # print("best result for {} is {}".format(k, sRes[0]))
        best_params[k] = sRes[0]
        # print("------------------------------")
    # print("------- best params ---------")
    # print(json.dumps(best_params, indent=2))
    # now update config with best params
    os.rename(glimmCFG, glimmCFG + ".init")
    with open(glimmCFG, "w") as outfile:
        for k, v in base_cfg.items():
            if k in best_params:
                outfile.write(f"{k} {best_params[k][0]}\n")
            else:
                outfile.write(f"{k} {v}\n")

    final_train = test_training(TrainData, test_models, tmpdir, tool="glimmerhmm")
    optimize_end = time.time()
    optimize_elapsed = time.strftime("%H:%M:%S", time.gmtime(optimize_end - train_end))
    log.info(
        "Initial training completed in {} and parameter optimization completed in {}s\n{}".format(
            train_elapsed,
            optimize_elapsed,
            json.dumps(final_train, indent=2, default=str),
        )
    )

    return {"location": TrainData, "train_results": final_train}


def run_glimmerhmm(genome, train_data, predictions, folder="/tmp", log=sys.stderr.write):
    """
    Run GlimmerHMM on the input genome to predict genes and output the results in GFF3 format.

    This function executes GlimmerHMM using the specified genome and training data, processes the raw output,
    and converts it to a GFF3 format file containing gene predictions.

    Args:
        genome (str): Path to the input genome file.
        train_data (str): Path to the training data file for GlimmerHMM.
        predictions (str): Path to save the predicted gene annotations in GFF3 format.
        folder (str, optional): Temporary folder to store intermediate files. Defaults to "/tmp".
        log (function, optional): Logger function to write debug and error messages. Defaults to sys.stderr.write.

    Returns:
        int: Returns 0 if successful.
    """

    def _glimmer2gff3(input, output):
        """
        scaffold_39     GlimmerHMM      mRNA    23692   25015   .       +       .       ID=scaffold_39.path1.gene12;Name=scaffold_39.path1.gene12
        scaffold_39     GlimmerHMM      CDS     23692   23886   .       +       0       ID=scaffold_39.cds12.1;Parent=scaffold_39.path1.gene12;Name=scaffold_39.path1.gene12;Note=initial-exon
        scaffold_39     GlimmerHMM      CDS     24282   24624   .       +       0       ID=scaffold_39.cds12.2;Parent=scaffold_39.path1.gene12;Name=scaffold_39.path1.gene12;Note=internal-exon
        scaffold_39     GlimmerHMM      CDS     24711   25015   .       +       2       ID=scaffold_39.cds12.3;Parent=scaffold_39.path1.gene12;Name=scaffold_39.path1.gene12;Note=final-exon
        scaffold_39     GlimmerHMM      mRNA    25874   27899   .       -       .       ID=scaffold_39.path1.gene13;Name=scaffold_39.path1.gene13
        scaffold_39     GlimmerHMM      CDS     25874   26973   .       -       2       ID=scaffold_39.cds13.1;Parent=scaffold_39.path1.gene13;Name=scaffold_39.path1.gene13;Note=final-exon
        scaffold_39     GlimmerHMM      CDS     27257   27899   .       -       0       ID=scaffold_39.cds13.2;Parent=scaffold_39.path1.gene13;Name=scaffold_39.path1.gene13;Note=initial-exon
        """
        with open(output, "w") as outfile:
            outfile.write("##gff-version 3\n")
            exonCounts = {}
            GeneCount = 1
            skipList = []
            idsSeen = {}
            with open(input, "r") as infile:
                for i, line in enumerate(infile):
                    if line.startswith("##sequence-region"):
                        idsSeen = {}
                    if line.startswith("#") or line.startswith("\n"):
                        continue
                    line = line.strip()
                    if line.count("\t") < 8:
                        log(
                            "ERROR parsing GlimmerHMM Raw output in line {}:\n   {}".format(
                                i + 1, line
                            )
                        )
                        continue
                    (
                        contig,
                        source,
                        feature,
                        start,
                        end,
                        score,
                        strand,
                        phase,
                        attributes,
                    ) = line.split("\t")
                    ID, Parent, Name = (None,) * 3
                    info = attributes.split(";")
                    for x in info:
                        if x.startswith("ID="):
                            ID = x.replace("ID=", "")
                        elif x.startswith("Parent="):
                            Parent = x.replace("Parent=", "")
                    if Parent and Parent in skipList:
                        continue
                    if feature == "mRNA":
                        genelen = int(end) - int(start)
                        if genelen < 150:
                            if ID not in skipList:
                                skipList.append(ID)
                            continue
                        geneID = f"{contig}-g{GeneCount}"
                        transID = f"{contig}-g{GeneCount}.t1"
                        idsSeen[ID] = (geneID, transID)
                        outfile.write(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tID={};Alias={};\n".format(
                                contig,
                                "glimmerhmm",
                                "gene",
                                start,
                                end,
                                score,
                                strand,
                                phase,
                                geneID,
                                ID,
                            )
                        )
                        outfile.write(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tID={};Parent={};Alias={};\n".format(
                                contig,
                                "glimmerhmm",
                                "mRNA",
                                start,
                                end,
                                ".",
                                strand,
                                ".",
                                transID,
                                geneID,
                                ID,
                            )
                        )
                        GeneCount += 1
                    elif feature == "CDS":
                        if Parent in idsSeen:
                            geneID, transID = idsSeen.get(Parent)
                            if transID not in exonCounts:
                                exonCounts[transID] = 1
                            else:
                                exonCounts[transID] += 1
                            num = exonCounts.get(transID)
                            outfile.write(
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tID={}.exon{};Parent={};\n".format(
                                    contig,
                                    "glimmerhmm",
                                    "exon",
                                    start,
                                    end,
                                    ".",
                                    strand,
                                    ".",
                                    transID,
                                    num,
                                    transID,
                                )
                            )
                            outfile.write(
                                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tID={}.cds;Parent={};\n".format(
                                    contig,
                                    "glimmerhmm",
                                    feature,
                                    start,
                                    end,
                                    score,
                                    strand,
                                    phase,
                                    transID,
                                    transID,
                                )
                            )
                        else:
                            log(
                                "ERROR parsing GlimmerHMM Raw output in line {}:\n   {}".format(
                                    i + 1, line
                                )
                            )

    glimmerRaw = os.path.join(
        folder,
        f"{os.path.basename(genome)}.{str(uuid.uuid4())[:8]}.glimmerHMM.output.raw",
    )
    # glimmer write a temp file so if running multiprocessing this is a problem
    # so run it from a unique temp location
    # also note that running glimmerhmm directly like this can only run a single contig
    cmd = [
        which_path("glimmerhmm"),
        os.path.abspath(genome),
        os.path.abspath(train_data),
        "-g",
        "-f",
        "-n",
        "1",
    ]
    with tempfile.TemporaryDirectory() as tmpdirname:
        runSubprocess(cmd, log, cwd=tmpdirname, stdout=glimmerRaw)
    # now convert to proper GFF3 format
    _glimmer2gff3(glimmerRaw, predictions)

    if os.path.isfile(glimmerRaw):
        os.remove(glimmerRaw)

    return 0


def train_genemark(
    genome,
    train_gff,
    test_models,
    fungus=False,
    folder="/tmp",
    cpus=1,
    log=sys.stderr.write,
):
    """
    Train a GeneMark model using a self-training approach.

    This function sets up a training directory, runs GeneMark training, and evaluates the model using test data.
    It supports fungal genomes and utilizes multiple CPUs for faster processing.

    Args:
        genome (str): Path to the genome file for training.
        train_gff (str): Path to the training GFF file.
        test_models (dict): Dictionary containing test models for evaluation.
        fungus (bool, optional): Indicates if the genome belongs to a fungus. Defaults to False.
        folder (str, optional): Directory path for storing temporary files. Defaults to "/tmp".
        cpus (int, optional): Number of CPUs to use for training. Defaults to 1.
        log (function, optional): Logger for writing debug and error messages. Defaults to sys.stderr.write.

    Returns:
        dict: A dictionary with the location of the trained GeneMark model and training results.
    """
    # self training with genemark
    # generate training directory ontop of the dir that is passed
    tmpdir = os.path.join(folder, "genemark")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    os.makedirs(tmpdir)

    # define training model
    genemark_mod = os.path.join(tmpdir, "gmhmm.mod")
    if os.path.abspath(genome) != os.path.abspath(os.path.join(folder, os.path.basename(genome))):
        shutil.copyfile(genome, os.path.join(folder, os.path.basename(genome)))
    cmd = [
        "gmes_petap.pl",
        "--ES",
        "--cores",
        str(cpus),
        "--work_dir",
        "genemark",
        "--training",
    ]
    if fungus:
        cmd.append("--fungus")
    cmd += ["--sequence", os.path.basename(genome)]
    train_start = time.time()
    runSubprocess(cmd, log, cwd=folder)
    if checkfile(genemark_mod):
        # then it appears to work, now we can test it
        init_train = test_training(
            genemark_mod,
            test_models,
            tmpdir,
            tool="genemark",
            log=log,
        )
        train_end = time.time()
        train_elapsed = time.strftime("%H:%M:%S", time.gmtime(train_end - train_start))
        log.info(
            "Initial training completed in {}s\n{}".format(
                train_elapsed, json.dumps(init_train, indent=2, default=str)
            )
        )
        return {"location": os.path.abspath(genemark_mod), "train_results": init_train}
    else:
        return {"location": None, "train_results": None}


def run_genemark(genome, train_data, predictions, folder="/tmp", log=sys.stderr.write):
    """
    Run GeneMark to predict genes in a genome using specified training data and save the predictions to a file.

    This function executes GeneMark on the provided genome and training data, processes the output,
    and converts it to GFF3 format for gene annotations.

    Args:
        genome (str): Path to the genome file.
        train_data (str): Path to the training data file.
        predictions (str): Path to save the predicted gene annotations.
        folder (str, optional): Directory to store temporary files. Defaults to "/tmp".
        log (function, optional): Logger function to write debug and error messages. Defaults to sys.stderr.write.

    Returns:
        None
    """
    # run genemark, note gmhmme3 only works on a single contig
    # and need to have both genome and training file in working directory
    tmpout = os.path.basename(predictions) + ".gtf"
    cmd = [
        "gmhmme3",
        "-f",
        "gtf",
        "-k",
        "0.03",
        "-m",
        os.path.basename(train_data),
        "-o",
        tmpout,
        os.path.basename(genome),
    ]
    workdir = os.path.dirname(genome)
    genemark_mod = os.path.join(workdir, os.path.basename(train_data))
    if not checkfile(genemark_mod):
        shutil.copy(os.path.abspath(train_data), genemark_mod)
    runSubprocess(cmd, log, cwd=workdir)
    Genes = gtf2dict(os.path.join(workdir, tmpout), genome)
    Cleaned = {}
    for k, v in natsorted(Genes.items()):
        new_gene = f"{v['contig']}-genemark.{k.replace('_g', '')}"
        new_model = v
        new_model["ids"] = [f"{new_gene}-T1"]
        Cleaned[new_gene] = new_model
    dict2gff3(Cleaned, output=predictions)


def train_snap(
    genome,
    train_gff,
    test_models,
    species=None,
    minintron=10,
    maxintron=3000,
    folder="/tmp",
    log=sys.stderr.write,
):
    """
    Train a SNAP model using the provided genome and gene models for training.

    This function prepares the training data, converts it to ZFF format, and runs SNAP training. It evaluates the model using test data and logs the results.

    Args:
        genome (str): Path to the genome file.
        train_gff (str or dict): Path to the GFF file containing gene models for training or a dictionary of gene models.
        test_models (dict): Dictionary of test gene models.
        species (str, optional): Species name. Defaults to None.
        minintron (int, optional): Minimum intron length. Defaults to 10.
        maxintron (int, optional): Maximum intron length. Defaults to 3000.
        folder (str, optional): Folder path for temporary files. Defaults to "/tmp".
        log (function, optional): Logger for debug and error messages. Defaults to sys.stderr.write.

    Returns:
        dict: A dictionary containing the location of the trained model and training results.
    """

    def _dict2zff(scaffoldDict, GeneDict, output):
        # take funannotate dictionary convert to zff training format
        with open(output, "w") as outfile:
            for k, v in natsorted(list(scaffoldDict.items())):
                outfile.write(">{:}\n".format(k))
                for genes in v:
                    gd = GeneDict.get(genes)
                    for i in range(0, len(gd["ids"])):
                        for num, c in enumerate(gd["CDS"][i]):
                            if gd["strand"] == "+":
                                start = c[0]
                                end = c[1]
                            else:
                                start = c[1]
                                end = c[0]
                            if num == 0:
                                outfile.write(
                                    "Einit\t{:}\t{:}\t{:}\n".format(start, end, gd["ids"][i])
                                )
                            elif num == len(gd["CDS"][i]) - 1:
                                outfile.write(
                                    "Eterm\t{:}\t{:}\t{:}\n".format(start, end, gd["ids"][i])
                                )
                            else:
                                outfile.write(
                                    "Exon\t{:}\t{:}\t{:}\n".format(start, end, gd["ids"][i])
                                )

    # generate training directory ontop of the dir that is passed
    tmpdir = os.path.join(folder, "snap")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    os.makedirs(tmpdir)

    # setup training files and check if exists
    snapHMM = os.path.join(tmpdir, "snap-trained.hmm")
    train_start = time.time()
    # load gff3 into dictionary
    if isinstance(train_gff, str):
        Genes = gff2dict(train_gff, genome)
    elif isinstance(train_gff, dict):
        Genes = train_gff
    scaff2genes = {}

    # load genome into dicrionary
    SeqRecords = fasta2dict(genome)

    # sort the dictionary
    def _sortDict(d):
        return (d[1]["contig"], d[1]["location"][0])

    sGenes = sorted(iter(Genes.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    scaff2genes = {}
    for k, v in list(sortedGenes.items()):
        if not v["contig"] in scaff2genes:
            scaff2genes[v["contig"]] = [k]
        else:
            scaff2genes[v["contig"]].append(k)

    # get only scaffolds that have gene models for training
    log.debug(
        "{:} gene models to train snap on {:} scaffolds".format(len(sGenes), len(scaff2genes))
    )
    trainingFasta = os.path.join(tmpdir, "snap-training.scaffolds.fasta")
    with open(trainingFasta, "w") as outfile:
        for title, seq in SeqRecords.items():
            if title in list(scaff2genes.keys()):
                outfile.write(">{:}\n{:}\n".format(title, softwrap(seq)))

    # convert to ZFF format
    origzff = os.path.join(tmpdir, "snap.training.zff")
    _dict2zff(scaff2genes, Genes, origzff)

    # now run SNAP training
    cmd = [
        "fathom",
        os.path.basename(origzff),
        os.path.basename(trainingFasta),
        "-categorize",
        "1000",
        "-min-intron",
        str(minintron),
        "-max-intron",
        str(maxintron),
    ]
    runSubprocess(cmd, log, cwd=tmpdir)

    cmd = ["fathom", "uni.ann", "uni.dna", "-export", "1000", "-plus"]
    runSubprocess(cmd, log, cwd=tmpdir)

    cmd = ["forge", "-min-intron", "10", "-lcmask", "export.ann", "export.dna"]
    runSubprocess(cmd, log, cwd=tmpdir)

    cmd = ["perl", which_path("hmm-assembler.pl"), "-r", "snap-trained", "."]
    runSubprocess(cmd, log, stdout=snapHMM, cwd=tmpdir)

    # now we can test the trained model
    train_elapsed = time.strftime("%H:%M:%S", time.gmtime(time.time() - train_start))
    init_train = test_training(snapHMM, test_models, tmpdir, tool="snap")
    log.info(
        "Initial training completed in {}s\n{}".format(
            train_elapsed, json.dumps(init_train, indent=2, default=str)
        )
    )

    # if this does not fail, then return some info in dict format
    return {"location": os.path.abspath(snapHMM), "train_results": init_train}


def run_snap(genome, train_data, predictions, folder="/tmp", log=sys.stderr.write):
    """
    Run SNAP prediction on a genome using training data and save the predictions in GFF3 format.

    This function executes SNAP on the provided genome and training data, processes the raw output,
    and converts it to GFF3 format for gene annotations.

    Args:
        genome (str): Path to the input genome file.
        train_data (str): Path to the training data file.
        predictions (str): Path to save the predicted annotations in GFF3 format.
        folder (str, optional): Directory to store temporary files. Defaults to '/tmp'.
        log (function, optional): Logger function for debug and error messages. Defaults to sys.stderr.write.

    Returns:
        int: Statistics about the predicted genes.
    """

    def _rawgff2gff3(input, fasta, output):
        """
        Scaffold_14	SNAP	Einit	804252	804255	7.048	-	.	Scaffold_14-snap.140
        Scaffold_14	SNAP	Exon	804163	804179	7.593	-	.	Scaffold_14-snap.140
        Scaffold_14	SNAP	Exon	803941	804099	-2.039	-	.	Scaffold_14-snap.140
        Scaffold_14	SNAP	Exon	803381	803862	-18.305	-	.	Scaffold_14-snap.140
        Scaffold_14	SNAP	Exon	802730	803259	5.042	-	.	Scaffold_14-snap.140
        Scaffold_14	SNAP	Eterm	802118	802671	-25.335	-	.	Scaffold_14-snap.140
        Scaffold_14	SNAP	Einit	807525	807529	7.047	-	.	Scaffold_14-snap.141
        Scaffold_14	SNAP	Exon	807420	807438	5.233	-	.	Scaffold_14-snap.141
        Scaffold_14	SNAP	Exon	807061	807342	49.005	-	.	Scaffold_14-snap.141
        Scaffold_14	SNAP	Exon	806749	807029	37.047	-	.	Scaffold_14-snap.141
        Scaffold_14	SNAP	Exon	805337	805357	11.064	-	.	Scaffold_14-snap.141
        Scaffold_14	SNAP	Eterm	805256	805307	19.709	-	.	Scaffold_14-snap.141
        Scaffold_14	SNAP	Einit	805508	805516	-1.215	+	.	Scaffold_14-snap.142
        Scaffold_14	SNAP	Exon	805537	806007	39.041	+	.	Scaffold_14-snap.142
        Scaffold_14	SNAP	Exon	806024	806389	20.790	+	.	Scaffold_14-snap.142
        Scaffold_14	SNAP	Exon	807668	807988	15.778	+	.	Scaffold_14-snap.142
        Scaffold_14	SNAP	Exon	811164	811182	11.201	+	.	Scaffold_14-snap.142
        Scaffold_14	SNAP	Exon	811230	811265	10.856	+	.	Scaffold_14-snap.142
        Scaffold_14	SNAP	Exon	812078	812083	16.816	+	.	Scaffold_14-snap.142
        Scaffold_14	SNAP	Exon	812322	812385	18.924	+	.	Scaffold_14-snap.142
        Scaffold_14	SNAP	Exon	813339	813414	11.856	+	.	Scaffold_14-snap.142
        Scaffold_14	SNAP	Exon	813886	813966	20.982	+	.	Scaffold_14-snap.142
        Scaffold_14	SNAP	Exon	814631	814799	27.101	+	.	Scaffold_14-snap.142
        Scaffold_14	SNAP	Eterm	814837	815456	35.426	+	.	Scaffold_14-snap.142
        """
        # need to load/generate a funannotate dictionary, then output to gff3 format
        Genes = {}
        contig = ""
        with open(input, "r") as infile:
            for line in infile:
                line = line.strip()
                if line.startswith("#") or line.startswith("\n"):
                    continue
                # basic gff easily parsable format here
                cols = line.split("\t")
                contig = cols[0]
                ID = cols[8]
                start = int(cols[3])
                end = int(cols[4])
                strand = cols[6]
                # score = cols[5]  # Unused variable
                phase = "?"
                if ID not in Genes:
                    Genes[ID] = {
                        "name": None,
                        "type": ["mRNA"],
                        "transcript": [],
                        "cds_transcript": [],
                        "protein": [],
                        "5UTR": [[]],
                        "3UTR": [[]],
                        "codon_start": [[]],
                        "ids": [ID + "-T1"],
                        "CDS": [[(start, end)]],
                        "mRNA": [[(start, end)]],
                        "strand": strand,
                        "location": (start, end),
                        "contig": contig,
                        "product": [[]],
                        "source": "snap",
                        "phase": [[phase]],
                        "db_xref": [[]],
                        "EC_number": [[]],
                        "gene_synonym": [],
                        "go_terms": [[]],
                        "note": [[]],
                        "partialStart": [[]],
                        "partialStop": [[]],
                        "pseudo": False,
                    }
                else:
                    Genes[ID]["CDS"][0].append((start, end))
                    Genes[ID]["mRNA"][0].append((start, end))
                    Genes[ID]["phase"][0].append(phase)
                    if start < Genes[ID]["location"][0]:
                        Genes[ID]["location"] = (start, Genes[ID]["location"][1])
                    if end > Genes[ID]["location"][1]:
                        Genes[ID]["location"] = (Genes[ID]["location"][0], end)
        # translate, check partial, etc
        SeqRecords = fasta2dict(fasta)
        # loop through and make sure CDS and exons are properly sorted and codon_start is correct, translate to protein space
        Genes = validate_models(Genes, SeqRecords, logger=log, table=1, gap_filter=False)
        # now write to GFF3
        dict2gff3(Genes, output)
        return 0

    # run snap here
    snapRaw = os.path.join(
        folder,
        f"{os.path.basename(genome)}.{str(uuid.uuid4())[:8]}.snap-prediction-raw.gff",
    )
    # now run SNAP prediction
    cmd = [
        "snap",
        "-gff",
        "-quiet",
        os.path.abspath(train_data),
        os.path.abspath(genome),
    ]
    with tempfile.TemporaryDirectory() as tmpdirname:
        runSubprocess(cmd, log, stdout=snapRaw, cwd=tmpdirname)
    # convert zff to proper gff3
    gene_stats = _rawgff2gff3(snapRaw, genome, predictions)
    # clean up raw
    if os.path.isfile(snapRaw):
        os.remove(snapRaw)
    return gene_stats


def run_augustus(
    genome,
    train_data,
    predictions,
    folder="/tmp",
    log=sys.stderr.write,
    hints=False,
    config_path=False,
):
    """
    Run Augustus gene prediction tool on the specified genome using the provided training data.

    This function executes Augustus with the specified parameters, processes the output, and saves the gene predictions in GFF3 format.

    Args:
        genome (str): Path to the genome file.
        train_data (str): Species-specific training data for Augustus.
        predictions (str): Path to store the predicted annotations in GFF3 format.
        folder (str, optional): Temporary folder to use for running Augustus. Defaults to '/tmp'.
        log (function, optional): Logger function for writing debug and error messages. Defaults to sys.stderr.write.
        hints (str, optional): Path to hints file for Augustus. Defaults to False.
        config_path (str, optional): Path to Augustus configuration file. Defaults to False.

    Returns:
        int: Returns 0 if successful.
    """
    # function to run augustus, can be slow if pass entire genome here
    cmd = [
        "augustus",
        f"--species={train_data}",
        "--gff3=on",
        "--UTR=off",
        "--stopCodonExcludedFromCDS=False",
        f"--extrinsicCfgFile={os.path.join(os.path.dirname(os.path.abspath(__file__)), 'resources', 'extrinsic.E.XNT.RM.cfg')}",
    ]
    if config_path:
        cmd.append(f"--AUGUSTUS_CONFIG_PATH={config_path}")
    if hints:
        cmd.append(f"--hintsfile={os.path.abspath(hints)}")
    cmd.append(os.path.abspath(genome))

    with tempfile.TemporaryDirectory() as tmpdirname:
        runSubprocess(cmd, log, stdout=predictions, cwd=tmpdirname)
    # finally parse the combined output and be done
    parse_augustus_gff3(predictions)

    return 0


def run_augustus_join(
    genome,
    train_data,
    predictions,
    folder="/tmp",
    log=sys.stderr.write,
    hints=False,
    config_path=False,
):
    """
    Run Augustus to predict genes on a genome using specified training data and hints, then join predictions.

    This function executes Augustus on the provided genome, processes the output in chunks to improve performance,
    joins the predictions, and parses the final output.

    Args:
        genome (str): Path to the genome file.
        train_data (str): Species-specific training data for Augustus.
        predictions (str): Path to store the predicted annotations in GFF3 format.
        folder (str, optional): Temporary folder to use for running Augustus. Defaults to '/tmp'.
        log (function, optional): Logger function for writing debug and error messages. Defaults to sys.stderr.write.
        hints (str, optional): Path to hints file for Augustus. Defaults to False.
        config_path (str, optional): Path to Augustus configuration file. Defaults to False.

    Returns:
        int: Returns 0 if successful.
    """
    # need to find the join_genes script
    JOINPREDS = which_path("join_aug_pred.pl")

    if not JOINPREDS:
        if env["AUGUSTUS_BASE"]:
            p = os.path.join(env["AUGUSTUS_BASE"], "join_aug_pred.pl")
            if os.path.isfile(p):
                JOINPREDS = p
    # function to run augustus, can be slow if pass entire genome here
    cmd = [
        "augustus",
        f"--species={train_data}",
        "--gff3=on",
        "--UTR=off",
        "--stopCodonExcludedFromCDS=False",
        f"--extrinsicCfgFile={os.path.join(os.path.dirname(os.path.abspath(__file__)), 'resources', 'extrinsic.E.XNT.RM.cfg')}",
    ]
    if config_path:
        cmd.append(f"--AUGUSTUS_CONFIG_PATH={config_path}")
    if hints:
        cmd.append(f"--hintsfile={os.path.abspath(hints)}")
    # get the prediction start/stop
    p_ranges = rangify_genome(genome)
    time_start = time.time()
    with tempfile.TemporaryDirectory() as tmpdirname:
        aug_preds = []
        for k, v in p_ranges.items():
            for i, r in enumerate(v):  # loop through each range as separate process
                aug_out = os.path.join(tmpdirname, f"{k}_part{i}.augustus.gff3")
                aug_preds.append(aug_out)
                run_cmd = cmd + [
                    f"--predictionStart={r[0]}",
                    f"--predictionEnd={r[1]}",
                    os.path.abspath(genome),
                ]
                runSubprocess(run_cmd, log, stdout=aug_out, cwd=tmpdirname)
        # when this is all finished, then need to join them all
        combined_preds = os.path.join(tmpdirname, "all_augustus_preds.gff3")
        with open(combined_preds, "w") as outfile:
            for fname in aug_preds:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
        # now join them all
        subprocess.run([JOINPREDS], stdin=combined_preds, stdout=predictions)
    # finally parse the combined output and be done
    parse_augustus_gff3(predictions)
    elapsed = time.time() - time_start
    elapsed_str = time.strftime("%Hh%Mm%Ss", time.gmtime(elapsed))
    if isinstance(log, types.BuiltinFunctionType):
        logger = log
    else:
        logger = log.debug
    logger(f"Augustus took {elapsed_str} to predict genes on {genome}")
    return 0


def rangify_genome(genome, overlap=20000):
    """
    Generate ranges for a given genome by breaking large contigs into chunks with optional overlap.

    This function processes the genome file, dividing large contigs into smaller chunks to improve processing speed,
    particularly for tools like Augustus that are slow with large contigs.

    Args:
        genome (str): Path to the genome file in FASTA format.
        overlap (int, optional): Overlap size between chunks. Defaults to 20000.

    Returns:
        dict: A dictionary containing ranges for each contig in the genome.
    """
    # augustus is really slow if given large contigs, break these up into chunks seems faster
    # this function will return dictionary of ranges
    ranges = {}
    for title, seq in pyfastx.Fasta(genome, build_index=False, uppercase=False):
        if len(seq) <= 1e6:
            ranges[title] = [(1, len(seq))]
        else:
            n_parts = int(len(seq) / 1e6 + 1)
            chunks = len(seq) / n_parts
            for i in range(0, n_parts):
                if i == 0:
                    start = 1
                    end = int(chunks + overlap)
                else:
                    start = int(end - overlap)
                    end = int(start + chunks + overlap)
                if end > len(seq):
                    end = len(seq)
                if title not in ranges:
                    ranges[title] = [(start, end)]
                else:
                    ranges[title].append((start, end))
    return ranges


def parse_augustus_gff3(aug, hiq_support=90):
    """
    Parse an Augustus GFF3 file to filter gene predictions with supporting evidence.

    This function processes an Augustus GFF3 file, identifies gene predictions with transcript support,
    and rewrites the GFF3 file with updated fields based on the level of support.

    Args:
        aug (str): Path to the Augustus GFF3 file.
        hiq_support (int, optional): Minimum percentage of transcript support for high-quality predictions. Defaults to 90.

    Returns:
        None
    """
    # this will look for gene predictions that have supporting evidence
    # will also re-write GFF3 properly
    bak = f"{aug}.bak"
    os.rename(aug, bak)
    with open(aug, "w") as outfile:
        outfile.write("##gff-version 3\n")
        with open(bak, "r") as augustus:
            for pred in readBlocks(augustus, "# start gene"):
                values = []
                geneID = ""
                support = ""
                gfflines = []
                if pred[0].startswith("# This output"):
                    continue
                if pred[0].startswith("##gff-version 3"):
                    continue
                for line in pred:
                    line = line.replace("\n", "")
                    if line.startswith("# start gene"):
                        geneID = line.split(" ")[-1]
                        values.append(geneID)
                    elif line.startswith("# % of transcript supported by hints"):
                        support = line.split(" ")[-1]
                        values.append(support)
                    elif not line.startswith("#"):
                        if "\tgene\t" in line or "\ttranscript\t" in line or "\tCDS\t" in line:
                            gfflines.append(line.split("\t"))
                ex_count = 0
                for x in gfflines:
                    x[5] = "."
                    x[8] = x[8].replace(geneID, f"{x[0]}-{geneID}")
                    if float(support) >= hiq_support:
                        x[1] = "augustus-hiq"
                    else:
                        x[1] = "augustus"
                    if x[2] == "transcript":
                        x[2] = "mRNA"
                    outfile.write("{}\n".format("\t".join(x)))
                    if x[2] == "CDS":
                        x[2] = "exon"
                        ex_count += 1
                        x[8] = x[8].replace(".cds", f".exon{ex_count}")
                        x[7] = "."
                        outfile.write("{}\n".format("\t".join(x)))
    os.remove(bak)


def train_augustus(
    genome,
    train_gff,
    test_models,
    species=None,
    minintron=10,
    maxintron=3000,
    folder="/tmp",
    log=sys.stderr.write,
    optimize=False,
    n_train=150,
    cpus=1,
):
    """
    Train an Augustus model using the provided genome and training data.

    Args:
        genome (str): Path to the genome file.
        train_gff (str): Path to the training GFF file.
        test_models (dict): Dictionary of test models.
        species (str, optional): Species name for training. Defaults to None.
        minintron (int, optional): Minimum intron length. Defaults to 10.
        maxintron (int, optional): Maximum intron length. Defaults to 3000.
        folder (str, optional): Folder path for temporary files. Defaults to "/tmp".
        log (function, optional): Logging function. Defaults to sys.stderr.write.
        optimize (bool or str, optional): Flag to optimize training. Defaults to False.
        n_train (int, optional): Number of training instances. Defaults to 150.
        cpus (int, optional): Number of CPUs to use. Defaults to 1.

    Returns:
        dict: Dictionary containing the location of the trained model, species name, and training results.
    """
    # get augustus training scripts
    OPTIMIZE = which_path("optimize_augustus.pl")
    NEW_SPECIES = which_path("new_species.pl")

    if not OPTIMIZE:
        if env["AUGUSTUS_BASE"]:
            p = os.path.join(env["AUGUSTUS_BASE"], "optimize_augustus.pl")
            if os.path.isfile(p):
                OPTIMIZE = p
    if not NEW_SPECIES:
        if env["AUGUSTUS_BASE"]:
            p = os.path.join(env["AUGUSTUS_BASE"], "new_species.pl")
            if os.path.isfile(p):
                NEW_SPECIES = p
    # need a check here if any of above are none
    for x in [OPTIMIZE, NEW_SPECIES]:
        if not x:
            raise Exception(
                "Missing Augustus training scripts, ensure that the augustus scripts directory is in your PATH"
            )

    if not species:
        species = str(uuid.uuid4())

    # isolate training files in separate folder
    tmpdir = os.path.join(folder, "augustus")
    if os.path.isdir(tmpdir):
        shutil.rmtree(tmpdir)
    os.makedirs(tmpdir)
    os.makedirs(os.path.join(tmpdir, "species"))

    train_start = time.time()
    # for augustus training needs to be in genbank format
    trainingset = os.path.join(tmpdir, "augustus.training.gb")

    # test gb_io
    gb_total, _ = gff3_to_gbio(
        genome, train_gff, trainingset, train_test=False, offset=1000, lowercase=True
    )

    # copy over augustus training files/dirs you'll need
    shutil.copytree(
        os.path.join(env["AUGUSTUS_CONFIG"], "species", "generic"),
        os.path.join(tmpdir, "species", "generic"),
    )
    shutil.copytree(
        os.path.join(env["AUGUSTUS_CONFIG"], "parameters"),
        os.path.join(tmpdir, "parameters"),
    )
    shutil.copytree(
        os.path.join(env["AUGUSTUS_CONFIG"], "model"),
        os.path.join(tmpdir, "model"),
    )
    shutil.copytree(
        os.path.join(env["AUGUSTUS_CONFIG"], "extrinsic"),
        os.path.join(tmpdir, "extrinsic"),
    )
    shutil.copytree(
        os.path.join(env["AUGUSTUS_CONFIG"], "profile"),
        os.path.join(tmpdir, "profile"),
    )
    # optimize from an existing training set
    train_log = os.path.join(tmpdir, "augustus_training.log")
    if optimize and optimize is not True:
        shutil.copytree(
            os.path.join(env["AUGUSTUS_CONFIG"], "species", optimize),
            os.path.join(tmpdir, "species", species),
        )
        # now need to rename the files
        files = []
        for f in os.listdir(os.path.join(tmpdir, "species", species)):
            files.append(os.path.join(tmpdir, "species", species, f))
        for x in files:
            if "parameters.cfg" in x:
                with open(x.replace(optimize, species), "w") as outfile:
                    with open(x, "r") as infile:
                        for line in infile:
                            outfile.write(line.replace(optimize, species))
                os.remove(x)
            else:
                os.rename(x, x.replace(optimize, species))
        # see what initial results are
        aug_init_training = os.path.join(tmpdir, "augustus.initial.training.txt")
        init_train = test_augustus_predictions(tmpdir, species, trainingset, aug_init_training)
        log.info(f"Augustus initial training results:\n{init_train}")
        # so idea here is to copy over existing and then run optimize augustus on it
        optimize_cmd = [
            OPTIMIZE,
            f"--AUGUSTUS_CONFIG_PATH={tmpdir}",
            f"--species={species}",
            "--UTR=off",
            f"--cpus={cpus}",
            "--rounds=1",
            f"{os.path.basename(trainingset)}",
        ]
        with open(train_log, "w") as logfile:
            log.debug(f"{' '.join(optimize_cmd)}")
            subprocess.run(optimize_cmd, cwd=tmpdir, stdout=logfile, stderr=logfile)
    else:
        new_species_cmd = [
            NEW_SPECIES,
            f"--AUGUSTUS_CONFIG_PATH={tmpdir}",
            f"--species={species}",
        ]

        etraining_cmd = [
            "etraining",
            f"--species={species}",
            f"--AUGUSTUS_CONFIG_PATH={tmpdir}",
            f"{os.path.basename(trainingset)}",
        ]
        with open(train_log, "w") as logfile:
            log.debug(f"{' '.join(new_species_cmd)}")
            subprocess.run(new_species_cmd, stdout=logfile, stderr=logfile)
            log.debug(f"{' '.join(etraining_cmd)}")
            subprocess.run(etraining_cmd, cwd=tmpdir, stdout=logfile, stderr=logfile)

        # see if optimize is True
        if optimize and optimize is True:
            aug_init_training = os.path.join(tmpdir, "augustus.initial.training.txt")
            init_train = test_augustus_predictions(tmpdir, species, trainingset, aug_init_training)
            log.info(f"Augustus initial training results:\n{init_train}")
            # so idea here is to copy over existing and then run optimize augustus on it
            optimize_cmd = [
                OPTIMIZE,
                f"--AUGUSTUS_CONFIG_PATH={tmpdir}",
                f"--species={species}",
                "--UTR=off",
                f"--cpus={cpus}",
                "--rounds=1",
                f"{os.path.basename(trainingset)}",
            ]
            with open(train_log, "a") as logfile:
                log.debug(f"{' '.join(optimize_cmd)}")
                subprocess.run(optimize_cmd, cwd=tmpdir, stdout=logfile, stderr=logfile)

    # now predict the test to see accuracy
    # now we can test the trained model
    init_train = test_training(species, test_models, tmpdir, tool="augustus", aug_config_dir=tmpdir)
    train_elapsed = time.strftime("%H:%M:%S", time.gmtime(time.time() - train_start))
    log.info(
        "Initial training completed in {}s\n{}".format(
            train_elapsed, json.dumps(init_train, indent=2, default=str)
        )
    )
    # aug_final_training = os.path.join(tmpdir, "augustus.final.training.txt")
    # final_train = test_augustus_predictions(
    #    tmpdir, species, trainingset, aug_final_training
    # )
    # log.info(f"Augustus final training results:\n{final_train}")

    # now return output
    return {
        "location": os.path.abspath(tmpdir),
        "species": species,
        "train_results": init_train,
    }


def test_augustus_predictions(config_dir, species, trainingset, output):
    """
    Test Augustus predictions against a reference set and output the comparison results.

    This function compares the predicted gene models from Augustus with a reference set of gene models,
    evaluates the accuracy, and writes the results to an output file.

    Args:
        predictions (str): Path to the file containing Augustus predictions in GFF3 format.
        reference (str): Path to the reference gene models file in GFF3 format.
        output (str): Path to save the comparison results.
        log (function, optional): Logger function for writing debug and error messages. Defaults to sys.stderr.write.

    Returns:
        None
    """
    aug_cmd = [
        "augustus",
        f"--AUGUSTUS_CONFIG_PATH={config_dir}",
        f"--species={species}",
        f"{trainingset}.test",
    ]
    with open(output, "w") as initrain:
        subprocess.run(aug_cmd, stdout=initrain)

    train_results = getTrainResults(output)
    trainTable = [
        ["Feature", "Specificity", "Sensitivity"],
        [
            "nucleotides",
            "{:.1%}".format(train_results[0]),
            "{:.1%}".format(train_results[1]),
        ],
        [
            "exons",
            "{:.1%}".format(train_results[2]),
            "{:.1%}".format(train_results[3]),
        ],
        [
            "genes",
            "{:.1%}".format(train_results[4]),
            "{:.1%}".format(train_results[5]),
        ],
    ]
    train_table = print_table(trainTable, return_str=True)
    return train_table


def getTrainResults(input):
    """
    Read a file containing training results and extract specific values related to nucleotide, exon, and gene levels.

    This function processes the input file to find lines starting with "nucleotide level", "exon level", and "gene level",
    extracts relevant values from these lines, and returns them as a tuple of floats.

    Args:
        input (str): Path to the file containing training results.

    Returns:
        tuple: A tuple containing six float values:
            - Two values from the nucleotide level line.
            - Two values from the exon level line.
            - Two values from the gene level line.
    """
    with open(input, "r") as train:
        for line in train:
            try:
                line = line.decode("utf-8")
            except AttributeError:
                pass
            line = line.rstrip()
            if line.startswith("nucleotide level"):
                line = line.replace(" ", "")
                values1 = line.split("|")  # get [1] and [2]
            if line.startswith("exon level"):
                line = line.replace(" ", "")  # get [6] and [7]
                values2 = line.split("|")
            if line.startswith("gene level"):
                line = line.replace(" ", "")
                values3 = line.split("|")  # get [6] and [7]
        return (
            float(values1[1]),
            float(values1[2]),
            float(values2[6]),
            float(values2[7]),
            float(values3[6]),
            float(values3[7]),
        )


def prot_alignments_to_hints(file):
    """
    Extract information from a GFF3 exonerate file to generate hints for gene prediction.

    This function processes the input file to create a contig-keyed dictionary with specific formatting rules,
    mimicking the exonerate2hints conversion. It adjusts CDS parts and introns to generate structured hints.

    Args:
        file (str): Path to the GFF3 exonerate file.

    Returns:
        dict: A dictionary where keys are contigs and values are lists of formatted hint lines.
    """
    # mimic exonerate2hints from GFF3 exonerate file
    # CDSpart +/- 15 bp to each match
    # intron as is
    # return a contig keyed dictionary with lines written as lists
    """
    #gff3 via EVM
    scaffold_20 exonerate   nucleotide_to_protein_match 225035  225823  82.13   +   .   ID=match.11677.2;Target=VC83_07547 1 96
    scaffold_20 exonerate   nucleotide_to_protein_match 53957   54342   92.93   +   .   ID=match.11677.3;Target=VC83_02595 1 129
    scaffold_20 exonerate   nucleotide_to_protein_match 54397   54904   92.93   +   .   ID=match.11677.3;Target=VC83_02595 130 299
    scaffold_107    exonerate   nucleotide_to_protein_match 77634   78119   89.95   -   .   ID=match.11677.5;Target=VC83_08471 1 163
    scaffold_107    exonerate   nucleotide_to_protein_match 77501   77546   89.95   -   .   ID=match.11677.5;Target=VC83_08471 163 178
    scaffold_107    exonerate   nucleotide_to_protein_match 77385   77422   89.95   -   .   ID=match.11677.5;Target=VC83_08471 179 191

    #corresponding exonerate2hints
    scaffold_20 xnt2h   CDSpart 225050  225808  .   +   .   src=XNT;grp=VC83_07547;pri=4
    scaffold_20 xnt2h   CDSpart 53972   54327   .   +   .   src=XNT;grp=VC83_02595;pri=4
    scaffold_20 xnt2h   intron  54343   54396   .   +   .   src=XNT;grp=VC83_02595;pri=4
    scaffold_20 xnt2h   CDSpart 54412   54889   .   +   .   src=XNT;grp=VC83_02595;pri=4
    scaffold_107    xnt2h   CDSpart 77649   78104   .   -   .   src=XNT;grp=VC83_08471;pri=4
    scaffold_107    xnt2h   intron  77547   77633   .   -   .   src=XNT;grp=VC83_08471;pri=4
    scaffold_107    xnt2h   CDSpart 77516   77531   .   -   .   src=XNT;grp=VC83_08471;pri=4
    scaffold_107    xnt2h   intron  77423   77500   .   -   .   src=XNT;grp=VC83_08471;pri=4
    scaffold_107    xnt2h   CDSpart 77400   77407   .   -   .   src=XNT;grp=VC83_08471;pri=4

    """
    Genes = {}
    with open(file, "r") as input:
        for line in input:
            if line.startswith("\n") or line.startswith("#"):
                continue
            line = line.rstrip()
            (
                contig,
                source,
                feature,
                start,
                end,
                score,
                strand,
                phase,
                attributes,
            ) = line.split("\t")
            start = int(start)
            end = int(end)
            ID, Target = (None,) * 2
            info = attributes.split(";")
            for x in info:
                if x.startswith("ID="):
                    ID = x.replace("ID=", "")
                elif x.startswith("Target="):
                    Target = x.replace("Target=", "").split(" ")[0]
            if ID not in Genes:
                Genes[ID] = {
                    "id": ID,
                    "target": Target,
                    "loc": [(start, end)],
                    "strand": strand,
                    "contig": contig,
                }
            else:
                Genes[ID]["loc"].append((start, end))
    # now lets sort through and write hints file
    Output = {}
    for k, v in natsorted(list(Genes.items())):
        if v["contig"] not in Output:
            Output[v["contig"]] = []
        if v["strand"] == "+":
            sortedCDS = sorted(v["loc"], key=lambda tup: tup[0])
            for i, x in enumerate(sortedCDS):  # loop through tuples
                Output[v["contig"]].append(
                    [
                        v["contig"],
                        "xnt2h",
                        "CDSpart",
                        x[0] - 15,
                        x[1] + 15,
                        ".",
                        v["strand"],
                        ".",
                        f"src=XNT;grp={v['target']};pri=4",
                    ]
                )
                if len(sortedCDS) > 1:
                    try:
                        Output[v["contig"]].append(
                            [
                                v["contig"],
                                "xnt2h",
                                "intron",
                                x[1] + 1,
                                sortedCDS[i + 1][0] - 1,
                                ".",
                                v["strand"],
                                ".",
                                f"src=XNT;grp={v['target']};pri=4",
                            ]
                        )
                    except IndexError:
                        pass
        else:
            sortedCDS = sorted(v["loc"], key=lambda tup: tup[0], reverse=True)
            for i, x in enumerate(sortedCDS):  # loop through tuples
                Output[v["contig"]].append(
                    [
                        v["contig"],
                        "xnt2h",
                        "CDSpart",
                        x[0] + 15,
                        x[1] - 15,
                        ".",
                        v["strand"],
                        ".",
                        f"src=XNT;grp={v['target']};pri=4",
                    ]
                )
                if len(sortedCDS) > 1:
                    try:
                        Output[v["contig"]].append(
                            [
                                v["contig"],
                                "xnt2h",
                                "intron",
                                sortedCDS[i + 1][1] + 1,
                                x[0] - 1,
                                ".",
                                v["strand"],
                                ".",
                                f"src=XNT;grp={v['target']};pri=4",
                            ]
                        )
                    except IndexError:
                        pass
    return Output


def alignments2dict(input, Genes):
    """
    Parse a transcript_alignments file and create a dictionary structure for each alignment.

    This function reads a transcript_alignments file, processes each line to extract relevant information,
    and populates a dictionary with alignment details such as mRNA coordinates, strand, percent identity,
    and additional attributes.

    Args:
        input (str): Path to the transcript_alignments file.
        Genes (dict): Dictionary to store the gene information.

    Returns:
        dict: Updated Genes dictionary with alignment information.
    """
    with open(input, "r") as infile:
        for line in infile:
            if line.startswith("\n") or line.startswith("#"):
                continue
            line = line.rstrip()
            (
                contig,
                source,
                feature,
                start,
                end,
                score,
                strand,
                phase,
                attributes,
            ) = line.split("\t")
            start = int(start)
            end = int(end)
            ID, Target, Extra = (None,) * 3
            for x in attributes.split(";"):
                if x.startswith("ID="):
                    ID = x.replace("ID=", "")
                elif x.startswith("Target="):
                    Target, Extra = x.split(" ", 1)
                    Target = Target.replace("Target=", "")
            if not ID:
                continue
            if ID not in Genes:
                Genes[ID] = {
                    "mRNA": [(start, end)],
                    "strand": strand,
                    "pident": [score],
                    "location": (start, end),
                    "contig": contig,
                    "extra": [Extra],
                }
            else:
                if contig != Genes[ID]["contig"]:
                    continue
                elif strand != Genes[ID]["strand"]:
                    continue
                else:
                    Genes[ID]["mRNA"].append((start, end))
                    Genes[ID]["pident"].append(score)
                    Genes[ID]["extra"].append(Extra)
                    # double check mRNA features are contained in gene coordinates
                    if start < Genes[ID]["location"][0]:
                        Genes[ID]["location"] = (start, Genes[ID]["location"][1])
                    if end > Genes[ID]["location"][1]:
                        Genes[ID]["location"] = (Genes[ID]["location"][0], end)
    return Genes


def introns_from_exons(input):
    """
    Generate a list of introns based on a list of exons.

    This function processes a list of exon tuples, where each tuple contains the start and end positions of an exon,
    and calculates the intron positions between consecutive exons.

    Args:
        input (list): A list of tuples representing exons, where each tuple contains the start and end positions of an exon.

    Returns:
        list: A list of tuples representing introns, where each tuple contains the start and end positions of an intron.
    """
    introns = []
    if len(input) > 1:
        for x, y in enumerate(input):
            try:
                introns.append((y[1] + 1, input[x + 1][0] - 1))
            except IndexError:
                pass
    return introns


def dict2hints(input):
    """
    Generate an Augustus hints file from a simple alignments dictionary.

    This function takes a dictionary containing alignment information, sorts the annotations by contig and start location,
    and generates hints for Augustus by processing exons and introns.

    Args:
        input (dict): A dictionary containing alignment information.

    Returns:
        dict: A dictionary representing the Augustus hints file, where keys are contigs and values are lists of hint lines.
    """
    from collections import OrderedDict

    def _sortDict(d):
        return (d[1]["contig"], d[1]["location"][0])

    # sort the annotations by contig and start location
    sGenes = natsorted(iter(input.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    Output = {}
    for k, v in list(sortedGenes.items()):
        if v["contig"] not in Output:
            Output[v["contig"]] = []
        sortedExons = sorted(v["mRNA"], key=lambda tup: tup[0])
        introns = introns_from_exons(sortedExons)
        for i, exon in enumerate(sortedExons):
            if i == 0 or i == len(sortedExons) - 1:
                Output[v["contig"]].append(
                    [
                        v["contig"],
                        "b2h",
                        "ep",
                        exon[0],
                        exon[1],
                        0,
                        v["strand"],
                        ".",
                        f"grp={k};pri=4;src=E",
                    ]
                )
            else:
                Output[v["contig"]].append(
                    [
                        v["contig"],
                        "b2h",
                        "exon",
                        exon[0],
                        exon[1],
                        0,
                        v["strand"],
                        ".",
                        f"grp={k};pri=4;src=E",
                    ]
                )
        if len(introns) > 0:
            for z in introns:
                Output[v["contig"]].append(
                    [
                        v["contig"],
                        "b2h",
                        "intron",
                        z[0],
                        z[1],
                        1,
                        v["strand"],
                        ".",
                        f"grp={k};pri=4;src=E",
                    ]
                )
    return Output


def evidence2hints(prot, tran, all_contigs, outdir):
    """
    Generate hints files for Augustus based on protein and transcript alignments for each contig.

    This function processes protein and transcript alignment files to create hints for gene prediction.
    It combines and sorts the hints for each contig and writes them to separate files, utilizing the
    Augustus script `join_mult_hints.pl` to merge the hints.

    Args:
        prot (str): Path to the protein alignments file.
        tran (str): Path to the transcript alignments file.
        all_contigs (list): List of paths to all contig FASTA files.
        outdir (str): Path to the output directory to store the hints files.

    Returns:
        None
    """
    # get augustus script
    JOINHINTS = which_path("join_mult_hints.pl")
    if not JOINHINTS:
        if env["AUGUSTUS_BASE"]:
            p = os.path.join(env["AUGUSTUS_BASE"], "join_mult_hints.pl")
            if os.path.isfile(p):
                JOINHINTS = p
    # parse the alignments, combine, sort, write hints to file per contig
    pHints = {}
    tHints = {}
    if checkfile(prot):
        pHints = prot_alignments_to_hints(prot)
    if checkfile(tran):
        t = {}
        t = alignments2dict(tran, t)
        tHints = dict2hints(t)
    for contig in all_contigs:
        contig = os.path.basename(contig).replace(".fasta", "")
        cHints = []
        if contig in pHints:
            cHints += pHints[contig]
        if contig in tHints:
            cHints += tHints[contig]
        # now sort and write, pipe output to augustus join_multi_hints.pl
        if len(cHints) > 0:
            finalhints = os.path.join(outdir, f"{contig}.fasta.hintsfile")
            with open(finalhints, "w") as outfile:
                p = subprocess.Popen([JOINHINTS], stdin=subprocess.PIPE, stdout=outfile)
                for h in natsorted(cHints, key=lambda x: (x[3], x[4], x[2], x[1])):
                    try:
                        p.stdin.write("{}\n".format("\t".join([str(x) for x in h])).encode())
                    except IOError as e:
                        if e.errno == errno.EPIPE or e.errno == errno.EINVAL:
                            # Stop loop on "Invalid pipe" or "Invalid argument".
                            # No sense in continuing with broken pipe.
                            break
                        else:
                            # Raise any other error.
                            raise
                p.stdin.close()
                p.wait()


def calculate_similarity(query, reference):
    """
    Calculate annotation edit distance (AED) and other metrics between two sets of coordinates.

    This function computes the AED, sensitivity, and precision at both nucleotide and exon levels
    by comparing query and reference coordinate sets. It evaluates true positives, false negatives,
    and false positives to derive these metrics.

    AED = 1 - (SN + SP / 2)
    SN = fraction of ref predicted
    SP = fraction query overlapping the ref

    And also return bases, exon, and gene level metrics
    Sensitivity = TP / (TP+FN)
    Precision =   TP / (TP+FP)

    TP = true positives: query bases/features which agree with ref bases/features
    FN = false negatives: reference bases/features which are missing in the query bases/features
    FP = false positives: query bases/features which are missing in the reference bases/features

    Args:
        query (list of tuples): List of tuples representing query coordinates.
        reference (list of tuples): List of tuples representing reference coordinates.

    Returns:
        tuple: A tuple containing the following metrics:
            - AED (float): Annotation Edit Distance.
            - nuc_sense (float): Nucleotide level sensitivity.
            - nuc_prec (float): Nucleotide level precision.
            - exon_sense (float): Exon level sensitivity.
            - exon_prec (float): Exon level precision.
    """

    def _length(listTup):
        total_len = 0
        for i in listTup:
            segment_len = abs(i[0] - i[1])
            total_len += segment_len
        return total_len

    # check if identical
    if query == reference:
        return 0.000, 1.00, 1.00, 1.00, 1.00
    # make sure sorted
    rLen = _length(reference)
    qLen = _length(query)
    refInterlap = InterLap(reference)
    query_exon_perf = []
    ref_exon_perf = []
    FPbases = 0
    FNbases = 0
    TPbases = 0
    refSeen = []
    refPerfect = []
    for exon in query:
        exon_perfect = False
        if exon in refInterlap:  # exon overlaps at least partially with reference
            hit = list(refInterlap.find(exon))
            for h in hit:
                refSeen.append(h)
                # will return array of exon minus hit at each pos
                diff = np.subtract(exon, h)
                if diff[0] == 0 and diff[1] == 0:  # then exon is perfect match
                    TPbases += abs(h[0] - h[1])
                    exon_perfect = True
                    refPerfect.append(h)
                elif diff[0] <= 0 and diff[1] >= 0:  # then query exon covers ref exon
                    cov = abs(h[0] - h[1])
                    FPbases += abs(diff[0])
                    FPbases += abs(diff[1])
                    TPbases += cov
                elif diff[0] <= 0 and diff[1] < 0:  # means query partial covers ref
                    cov = abs(h[0] - exon[1])
                    FPbases += abs(diff[0])
                    FNbases += abs(diff[1])
                    TPbases += cov
                elif diff[0] > 0 and diff[1] >= 0:  # means query partial covers ref
                    cov = abs(exon[0] - h[1])
                    FPbases += abs(diff[1])
                    FNbases += abs(diff[0])
                    TPbases += cov
                elif diff[0] > 0 and diff[1] < 1:
                    cov = abs(exon[0] - exon[1])
                    FNbases += abs(diff[1])
                    FNbases += abs(diff[0])
                    TPbases += cov
        else:
            FPbases += abs(exon[0] - exon[1])
        query_exon_perf.append(exon_perfect)
    # last thing is to double check reference for missed exons
    for ex in refInterlap:
        rex_perf = False
        if ex not in refSeen:
            FNbases += abs(ex[0] - ex[1])
        if ex in refPerfect:
            rex_perf = True
        ref_exon_perf.append(rex_perf)
    # calculate stats
    nuc_sense = TPbases / (TPbases + FNbases)
    nuc_prec = TPbases / (TPbases + FPbases)
    # calculate exons stats
    eTP = 0
    eFP = 0
    eFN = 0
    for x in query_exon_perf:
        if x:
            eTP += 1
        else:
            eFP += 1
    for y in ref_exon_perf:
        if y:
            eFN += 1
    # calculate exon stats
    try:
        exon_sense = eTP / (eTP + eFN)
    except ZeroDivisionError:
        exon_sense = 0.00
    try:
        exon_prec = eTP / (eTP + eFP)
    except ZeroDivisionError:
        exon_prec = 0.00
    # calc AED
    SN = (rLen - FNbases) / rLen
    SP = (qLen - FPbases) / qLen
    try:
        AED = 1 - ((SN + SP) / 2)
    except ZeroDivisionError:
        AED = 1.00
    return AED, nuc_sense, nuc_prec, exon_sense, exon_prec


def test_training(
    train_data,
    test_models,
    tmpdir,
    tool="glimmerhmm",
    aug_config_dir=".",
    log=sys.stderr.write,
):
    """
    Test the model with training data that it has never seen before.

    This function evaluates the performance of a gene prediction model using unseen training data. It writes the model test DNA, runs the test on the test contigs using the specified tool, and calculates various metrics such as Annotation Edit Distance (AED), sensitivity, and precision at nucleotide and exon levels.

    Args:
        train_data (str): Path to the training data file.
        test_models (dict): Dictionary containing test models with their respective training DNA and coordinates.
        tmpdir (str): Path to the temporary directory for storing intermediate files.
        tool (str, optional): The gene prediction tool to use. Options are "glimmerhmm", "snap", "augustus", or "genemark". Defaults to "glimmerhmm".
        aug_config_dir (str, optional): Path to the Augustus configuration directory. Defaults to ".".
        log (function, optional): Logger function for writing debug and error messages. Defaults to sys.stderr.write.

    Returns:
        dict: A dictionary containing the evaluation results, including AED, nucleotide sensitivity, nucleotide precision, exon sensitivity, exon precision, gene sensitivity, and gene precision.
    """
    # goal here is to test the model with training data (its never seen before)
    # now we want to test the training, first write the model test dna out
    if not os.path.isdir(os.path.join(tmpdir, "test")):
        os.makedirs(os.path.join(tmpdir, "test"))
    contigs = []
    for k, v in test_models.items():
        dest = os.path.join(tmpdir, "test", f"{k}.fa")
        contigs.append(dest)
        if not os.path.isfile(dest):
            with open(dest, "w") as outfile:
                outfile.write(f">{k}\n{softwrap(v['train_dna'])}\n")

    # now run the test on the test contigs, note that contig is the gene model name
    preds = {}
    for c in contigs:
        if tool == "glimmerhmm":
            cmd = [
                which_path("glimmerhmm"),
                os.path.abspath(c),
                os.path.abspath(train_data),
                "-g",
                "-f",
                "-n",
                "1",
            ]
            with tempfile.TemporaryDirectory() as tmpdirname:
                for line in execute(cmd, cwd=tmpdirname):
                    if line.count("\t") == 8:
                        cols = line.rstrip().split("\t")
                        if cols[2] == "CDS":
                            start = int(cols[3])
                            end = int(cols[4])
                            strand = cols[6]
                            gID = None
                            for x in cols[8].split(";"):
                                if x.startswith("Parent="):
                                    gID = x.split("=")[-1]
                            if gID:
                                if cols[0] not in preds:
                                    preds[cols[0]] = {
                                        gID: {
                                            "coords": [(start, end)],
                                            "strand": strand,
                                        }
                                    }
                                else:
                                    if gID not in preds[cols[0]]:
                                        preds[cols[0]][gID] = {
                                            "coords": [(start, end)],
                                            "strand": strand,
                                        }
                                    else:
                                        preds[cols[0]][gID]["coords"].append((start, end))
        elif tool == "snap":
            cmd = ["snap", os.path.abspath(train_data), os.path.abspath(c)]
            with tempfile.TemporaryDirectory() as tmpdirname:
                for line in execute(cmd, cwd=tmpdirname):
                    line = line.strip()
                    if line.startswith("#") or line.startswith("\n"):
                        continue
                    elif line.startswith(">"):
                        contig = line[1:]
                    else:
                        (
                            feature,
                            start,
                            end,
                            strand,
                            score,
                            fiveo,
                            threeo,
                            phase,
                            ID,
                        ) = line.split("\t")
                        start = int(start)
                        end = int(end)
                        # phase = "?"  # Unused variable
                        if contig not in preds:
                            preds[contig] = {
                                ID: {
                                    "coords": [(start, end)],
                                    "strand": strand,
                                }
                            }
                        else:
                            if ID not in preds[contig]:
                                preds[contig][ID] = {
                                    "coords": [(start, end)],
                                    "strand": strand,
                                }
                            else:
                                preds[contig][ID]["coords"].append((start, end))
        elif tool == "augustus":
            cmd = [
                "augustus",
                f"--AUGUSTUS_CONFIG_PATH={os.path.abspath(aug_config_dir)}",
                f"--species={train_data}",
                "--gff3=on",
                "--UTR=off",
                "--stopCodonExcludedFromCDS=False",
                f"--extrinsicCfgFile={os.path.join(os.path.dirname(os.path.abspath(__file__)), 'resources', 'extrinsic.E.XNT.RM.cfg')}",
                os.path.abspath(c),
            ]
            with tempfile.TemporaryDirectory() as tmpdirname:
                for line in execute(cmd, cwd=tmpdirname):
                    if line.startswith("#"):
                        continue
                    if line.count("\t") == 8:
                        cols = line.rstrip().split("\t")
                        if cols[2] == "CDS":
                            start = int(cols[3])
                            end = int(cols[4])
                            strand = cols[6]
                            gID = None
                            for x in cols[8].split(";"):
                                if x.startswith("Parent="):
                                    gID = x.split("=")[-1]
                            if gID:
                                if cols[0] not in preds:
                                    preds[cols[0]] = {
                                        gID: {
                                            "coords": [(start, end)],
                                            "strand": strand,
                                        }
                                    }
                                else:
                                    if gID not in preds[cols[0]]:
                                        preds[cols[0]][gID] = {
                                            "coords": [(start, end)],
                                            "strand": strand,
                                        }
                                    else:
                                        preds[cols[0]][gID]["coords"].append((start, end))
        elif tool == "genemark":
            # need to copy over the mod file as has to be in same directory
            shutil.copyfile(train_data, os.path.join(tmpdir, "test", os.path.basename(train_data)))
            g_out = os.path.join(tmpdir, "test", os.path.basename(c) + ".gtf")
            cmd = [
                "gmhmme3",
                "-f",
                "gtf",
                "-k",
                "0.03",
                "-m",
                os.path.basename(train_data),
                "-o",
                os.path.basename(c) + ".gtf",
                os.path.basename(c),
            ]
            runSubprocess(cmd, log, cwd=os.path.join(tmpdir, "test"))
            if checkfile(g_out):
                with open(g_out) as infile:
                    for line in infile:
                        if line.startswith(("#", "\n")):
                            continue
                        if line.count("\t") == 8:
                            cols = line.rstrip().split("\t")
                            if cols[2] == "CDS":
                                start = int(cols[3])
                                end = int(cols[4])
                                strand = cols[6]
                                gID = None
                                for x in cols[8].split(";"):
                                    x = x.strip()
                                    if x.startswith("gene_id"):
                                        gID = x.split()[-1].replace('"', "")
                                if gID:
                                    if cols[0] not in preds:
                                        preds[cols[0]] = {
                                            gID: {
                                                "coords": [(start, end)],
                                                "strand": strand,
                                            }
                                        }
                                    else:
                                        if gID not in preds[cols[0]]:
                                            preds[cols[0]][gID] = {
                                                "coords": [(start, end)],
                                                "strand": strand,
                                            }
                                        else:
                                            preds[cols[0]][gID]["coords"].append((start, end))

    # preds is keyed by gene name, and the value is keyed by gene model name, and the value is a dict of coords and strand
    data = {
        "aed": [],
        "nucleotide_sensitivity": [],
        "nucleotide_precision": [],
        "exon_sensitivity": [],
        "exon_precision": [],
        "rogue": [],
    }
    for k, v in preds.items():
        expected_coords = test_models[k]["train_coords"]
        expected_strand = test_models[k]["strand"]
        if len(v) == 1:  # only 1 model predicted
            gID = list(v.keys())[0]
            if v[gID]["strand"] != expected_strand:
                data["rogue"].append(gID)
                continue
            (
                aed,
                nucl_sens,
                nucl_precise,
                exon_sens,
                exon_precise,
            ) = calculate_similarity(
                sorted(v[gID]["coords"], key=lambda x: x[0]),
                sorted(expected_coords, key=lambda x: x[0]),
            )
            data["aed"].append(aed)
            data["nucleotide_sensitivity"].append(nucl_sens)
            data["nucleotide_precision"].append(nucl_precise)
            data["exon_sensitivity"].append(exon_sens)
            data["exon_precision"].append(exon_precise)
        else:  # means more than 1 model, so now need to pick the best one
            best_aed = 1.00
            best_model = None
            best_data = {}
            for m, cds in v.items():
                if cds["strand"] != expected_strand:
                    data["rogue"].append(m)
                    continue
                (
                    aed,
                    nucl_sens,
                    nucl_precise,
                    exon_sens,
                    exon_precise,
                ) = calculate_similarity(
                    sorted(cds["coords"], key=lambda x: x[0]),
                    sorted(expected_coords, key=lambda x: x[0]),
                )
                best_data[m] = [aed, nucl_sens, nucl_precise, exon_sens, exon_precise]
                if aed < best_aed:
                    best_aed = aed
                    best_model = m
            for y, z in best_data.items():
                if y != best_model:
                    if y not in data["rogue"]:
                        data["rogue"].append(y)
                else:
                    data["aed"].append(z[0])
                    data["nucleotide_sensitivity"].append(z[1])
                    data["nucleotide_precision"].append(z[2])
                    data["exon_sensitivity"].append(z[3])
                    data["exon_precision"].append(z[4])
    # calculate gene level metrics
    gTP = len([x for x in data["aed"] if x == 0.00])
    gFN = len([x for x in list(test_models.keys()) if x not in preds])
    gFP = len(data["rogue"])
    try:
        gene_sens = gTP / (gTP + gFN)
    except ZeroDivisionError:
        gene_sens = 0.00
    try:
        gene_prec = gTP / (gTP + gFP)
    except ZeroDivisionError:
        gene_prec = 0.00
    # tidy up results in a dictionary
    result = {
        "tool": tool,
        "model": os.path.basename(train_data),
        "n_test_genes": len(preds) + gFN,
        "ref_genes_found": len(preds),
        "ref_genes_missed": gFN,
        "extra_query_genes": len(data["rogue"]),
        "average_aed": sum(data["aed"]) / len(data["aed"]),
        "nucleotide_sensitivity": sum(data["nucleotide_sensitivity"])
        / len(data["nucleotide_sensitivity"]),
        "nucleotide_precision": sum(data["nucleotide_precision"])
        / len(data["nucleotide_precision"]),
        "exon_sensitivity": sum(data["exon_sensitivity"]) / len(data["exon_sensitivity"]),
        "exon_precision": sum(data["exon_precision"]) / len(data["exon_precision"]),
        "gene_sensitivity": gene_sens,
        "gene_precision": gene_prec,
    }
    return result
