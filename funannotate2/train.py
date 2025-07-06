import datetime
import getpass
import json
import os
import random
import shutil
import sys
from collections import OrderedDict, defaultdict

from buscolite.busco import runbusco
from buscolite.gff import gffwriter
from buscolite.utilities import summary_writer
from gfftk.fasta import fasta2dict, softwrap
from gfftk.gff import dict2gff3, gff2dict
from natsort import natsorted

from .abinitio import train_augustus, train_genemark, train_glimmerhmm, train_snap
from .config import env
from .fastx import analyzeAssemblySimple, simplify_headers
from .interlap import InterLap
from .log import finishLogging, startLogging, system_info
from .utilities import (
    checkfile,
    choose_best_augustus_species,
    choose_best_busco_species,
    create_directories,
    download,
    load_json,
    lookup_taxonomy,
    naming_slug,
    runSubprocess,
    which_path,
    get_odb_version,
    rename_gff_contigs,
    validate_busco_lineage,
    validate_augustus_species,
)
from .config import busco_taxonomy, augustus_species


def train(args):
    """
    Train ab initio gene prediction tools using a provided training set or generate one using BUSCO analysis.

    This function sets up the environment for training gene prediction tools such as Augustus, SNAP, and GlimmerHMM.
    It begins by creating necessary directories and setting up logging. The genome is loaded, and quality checks are performed.
    If a training set is not provided, BUSCO analysis is run to generate one. The training set is filtered and split into
    test and train sets. The function then executes training routines for each tool and saves the results in a JSON file.

    Parameters:
    - args (argparse.Namespace): Command-line arguments passed to the function.

    Returns:
    None
    """
    # This function is a wrapper for training augustus, snap, glimmerhmm, etc. If a training set is not passed
    # then it will run buscolite to generate a training set. Training set in GFF3 format will then be used to run
    # the training functions located in abinitio.py, ie train_snap(), train_augustus(), train_glimmerhmm(), etc.
    misc_dir, res_dir, log_dir = create_directories(args.out, base="train")
    # start logger
    global logger
    logger = startLogging(logfile=os.path.join(log_dir, "funannotate-train.log"))
    log = logger.info
    system_info(log)

    # load genome and do some QC checks
    # for downstream processing lets output the original genome to train_results
    original_genome = os.path.join(res_dir, os.path.basename(args.fasta))
    shutil.copyfile(args.fasta, original_genome)
    # set max header high as going to simplify below if no other issues
    logger.info("Loading genome assembly and running QC checks")
    stats, bad_names, nuc_errors = analyzeAssemblySimple(
        args.fasta, header_max=args.header_length
    )
    if len(bad_names) > 0:
        bad_string = ", ".join(bad_names)
        logger.critical(
            f"{len(bad_names)} contigs have headers that are longer than maximum of {args.header_length} characters\n{bad_string}"
        )
        raise SystemExit(1)
    if len(nuc_errors) > 0:
        bad_string = ""
        for x in nuc_errors:
            bad_string += f"{x[0]}: {x[1]}, "
        logger.critical(
            f"{len(nuc_errors)} contigs contain non-IUPAC characters\n{bad_string}"
        )
        raise SystemExit(1)
    logger.info(f"Genome stats:\n{json.dumps(stats, indent=2)}")

    # to circumvent any downstream issues, rename headers
    GenomeFasta = os.path.join(misc_dir, "genome.fasta")
    contigHeaderMap = simplify_headers(args.fasta, GenomeFasta, base="contig_")

    # get taxonomy information
    taxonomy = lookup_taxonomy(args.species)
    logger.info(f"Getting taxonomy information\n{json.dumps(taxonomy, indent=2)}")

    # validate and set augustus species
    if args.augustus_species:
        if not validate_augustus_species(args.augustus_species):
            logger.critical(f"Invalid Augustus species: {args.augustus_species}")
            logger.critical(
                f"Valid options are: {', '.join(sorted(augustus_species.keys()))}"
            )
            raise SystemExit(1)
        aug_species = args.augustus_species
        logger.info(f"Using user-specified Augustus species: {aug_species}")
    else:
        # choose best augustus species based on taxonomy
        aug_species = choose_best_augustus_species(taxonomy)
        logger.info(f"Choosing best augustus species based on taxonomy: {aug_species}")

    # validate and set busco lineage
    if args.busco_lineage:
        if not validate_busco_lineage(args.busco_lineage):
            logger.critical(f"Invalid BUSCO lineage: {args.busco_lineage}")
            logger.critical(
                f"Valid options are: {', '.join(sorted(busco_taxonomy.keys()))}"
            )
            raise SystemExit(1)
        busco_species = args.busco_lineage
        logger.info(f"Using user-specified BUSCO lineage: {busco_species}")
    else:
        # choose best busco species
        busco_species = choose_best_busco_species(taxonomy)
    # pull the latest odb version from downloads link
    odb_version = get_odb_version(
        os.path.join(os.path.dirname(__file__), "downloads.json")
    )
    busco_model_path = os.path.join(
        env["FUNANNOTATE2_DB"], f"{busco_species}_{odb_version}"
    )

    # run buscolite on genome to get training set
    filt_train_models = os.path.join(misc_dir, "training-models.final.gff3")
    if not os.path.isfile(filt_train_models):
        if args.training_set is None:
            args.training_set = os.path.join(misc_dir, "busco_training_set.gff3")
            if not os.path.isfile(args.training_set):
                logger.info(
                    f"Choosing best busco species based on taxonomy: {busco_species}"
                )
                if not os.path.isdir(busco_model_path):
                    download_urls = load_json(
                        os.path.join(os.path.dirname(__file__), "downloads.json")
                    )
                    busco_url = download_urls["busco"][busco_species][0]
                    busco_tgz = os.path.join(
                        env["FUNANNOTATE2_DB"], os.path.basename(busco_url)
                    )
                    logger.info(
                        f"Downloading {busco_species}_{odb_version} model from {busco_url}"
                    )
                    download(busco_url, busco_tgz, wget=False)
                    if os.path.isfile(busco_tgz):
                        runSubprocess(
                            ["tar", "-zxf", os.path.basename(busco_tgz)],
                            logger,
                            cwd=env["FUNANNOTATE2_DB"],
                        )
                        if os.path.isdir(busco_model_path):
                            os.remove(busco_tgz)
                log("Running buscolite to generate training set")
                buscolite(
                    GenomeFasta,
                    busco_model_path,
                    args.training_set,
                    species=aug_species,
                    cpus=args.cpus,
                    log=logger,
                )
                if not checkfile(args.training_set):
                    logger.critical("BUSCOlite failed to generate training set")
                    raise SystemExit(1)
            else:
                logger.info(f"Existing BUSCO results found: {args.training_set}")
        else:
            # however, first we need to change the contig names to the temporary ones
            train_temp = os.path.join(misc_dir, "training-set.temp.gff3")
            # need to invert contigHeaderMap
            gffRenameMap = {v: k for k, v in contigHeaderMap.items()}
            rename_gff_contigs(args.training_set, train_temp, gffRenameMap)
            args.training_set = train_temp

        # load  GFF3 training set, load with gfftk
        raw_train_set = gff2dict(args.training_set, GenomeFasta)
        # user might pass in a training set that has other models, we only want protein coding, so filter
        train_set = {
            k: v
            for k, v in raw_train_set.items()
            if "mRNA" in v["type"] and len(v["CDS"][0]) > 0
        }
        logger.info(
            f"Training set [{args.training_set}] loaded with {len(train_set)} gene models"
        )
        if len(train_set) == 0:
            logger.critical(
                f"No gene models found in training set: {args.training_set}"
            )
            raise SystemExit(1)

        # now filter training set
        models4training = selectTrainingModels(GenomeFasta, train_set, tmpdir=misc_dir)
        dict2gff3(models4training, filt_train_models)
    else:
        logger.info(f"Using existing training set: {filt_train_models}")
        models4training = trainmodels2dict(GenomeFasta, filt_train_models)

    # split into test/train sets
    n_test = int(len(models4training) * 0.20)
    if n_test > 200:
        n_test = 200
    logger.info(
        f"{len(models4training)} gene models selected for training, now splitting into test [n={n_test}] and train [n={len(models4training) - n_test}]"
    )

    test_model_keys = random.sample(list(models4training.keys()), n_test)
    train_model_keys = [x for x in models4training.keys() if x not in test_model_keys]
    test_models = {k: models4training[k] for k in test_model_keys}
    train_models = {k: models4training[k] for k in train_model_keys}
    filt_train_models_final = os.path.join(misc_dir, "training-models.train.gff3")
    dict2gff3(train_models, filt_train_models_final)

    # run augustus training functions
    logger.info("Training augustus using training set")
    augustus_train = train_augustus(
        GenomeFasta,
        train_models,
        test_models,
        folder=misc_dir,
        log=logger,
        cpus=args.cpus,
        optimize=False,
    )
    augustus_train["training_set"] = filt_train_models_final

    # run snap training functions
    logger.info("Training snap using training set")
    snap_train = train_snap(
        GenomeFasta, train_models, test_models, folder=misc_dir, log=logger
    )
    snap_train["training_set"] = filt_train_models_final

    # run glimmerHMM training functions
    logger.info("Training glimmerHMM using training set")
    glimm_train = train_glimmerhmm(
        GenomeFasta,
        train_models,
        test_models,
        folder=misc_dir,
        log=logger,
        asm_size=stats["size"],
    )
    glimm_train["training_set"] = filt_train_models_final

    # if adding more software make sure to include here
    train_data = {
        "augustus": augustus_train,
        "glimmerhmm": glimm_train,
        "snap": snap_train,
    }

    # check if genemark is installed, if so then run self training
    if which_path("gmes_petap.pl") and which_path("gmhmme3"):
        logger.info("Training GeneMark-ES using self-training")
        fungus_flag = False
        if taxonomy["kingdom"] == "Fungi":
            fungus_flag = True
        genemark_train = train_genemark(
            GenomeFasta,
            train_models,
            test_models,
            folder=misc_dir,
            fungus=fungus_flag,
            cpus=args.cpus,
            log=logger,
        )
        genemark_train["training_set"] = "self training"
        train_data["genemark"] = genemark_train
    else:
        logger.info(
            f"GeneMark-ES not installed, skipping training: gmes_petap.pl in PATH={which_path('gmes_petap.pl')}; gmhmme3 in PATH={which_path('gmhmme3')}"
        )

    # now lets get these data together and save as parameters.json file
    name_slug = naming_slug(args.species, args.strain)
    params_json = os.path.join(res_dir, f"{name_slug}.params.json")
    tdata = {
        "name": name_slug,
        "species": args.species,
        "taxonomy": taxonomy,
        "busco-lineage": busco_model_path,
        "assembly": stats,
        "abinitio": train_data,
        "date": datetime.datetime.now(),
        "user": getpass.getuser(),
    }
    with open(params_json, "w") as outfile:
        outfile.write(json.dumps(tdata, indent=2, default=str))
    logger.info(f"Ab initio training finished: {params_json}")
    logger.info(
        "The params.json file can be passed to funannotate2 predict or installed globally with funannotate2 species"
    )

    # finish
    finishLogging(log, vars(sys.modules[__name__])["__name__"])


def buscolite(
    genome,
    lineage,
    output,
    species="anidulans",
    cpus=1,
    flanks=2000,
    mode="genome",
    log=sys.stderr.write,
    summary=False,
):
    """
    Run BUSCO analysis on a genome using the specified lineage.

    This function performs a BUSCO analysis on the provided genome file using the specified
    lineage. It writes the results to an output file and optionally generates a summary file.

    Parameters:
    - genome (str): Path to the input genome file.
    - lineage (str): BUSCO lineage to use for analysis.
    - output (str): Path to the output file to write the results.
    - species (str, optional): Species to use for analysis (default is "anidulans").
    - cpus (int, optional): Number of CPUs to use for analysis (default is 1).
    - flanks (int, optional): Offset value for flanking regions (default is 2000).
    - mode (str, optional): Analysis mode (default is "genome").
    - log (function, optional): Logger function (default is sys.stderr.write).
    - summary (bool, optional): Whether to generate a summary file (default is False).

    Returns:
    None
    """
    d, m, stats, cfg = runbusco(
        genome,
        lineage,
        species=species,
        mode=mode,
        cpus=cpus,
        offset=flanks,
        logger=log,
        check_augustus=False,
        verbosity=1,
    )
    log.info(
        "Analysis complete:\n single-copy={}\n fragmented={}\n duplicated={}\n total={}".format(
            stats["single-copy"],
            stats["fragmented"],
            stats["duplicated"],
            stats["total"],
        )
    )
    with open(output, "w") as outfile:
        gffwriter(d, outfile)
    if summary:
        with open(summary, "w") as outfile:
            summary_writer(d, m, [], cfg, mode=mode)


def getTrainResults(infile):
    """
    Parse training results from the input file and extract specific metrics.

    This function reads a file containing training results and extracts specific values
    from lines starting with "nucleotide level", "exon level", and "gene level". It returns
    these values as a tuple of floats.

    Parameters:
    - infile (str): Path to the input file containing training results.

    Returns:
    - tuple: A tuple containing the extracted float values:
      - Nucleotide level value at index 1
      - Nucleotide level value at index 2
      - Exon level value at index 6
      - Exon level value at index 7
      - Gene level value at index 6
      - Gene level value at index 7
    """
    with open(infile, "r") as train:
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


def count_multi_CDS_genes(indict):
    """
    Count genes with more than one CDS in a funannotate annotation dictionary.

    This function iterates over a dictionary of funannotate annotations and counts the number
    of genes that have more than one coding sequence (CDS). It returns the total number of genes
    and the number of genes with multiple CDS.

    Parameters:
    - indict (dict): A dictionary containing funannotate annotations.

    Returns:
    - tuple: A tuple containing:
      - int: Total number of genes.
      - int: Number of genes with more than one CDS.
    """
    # take funannotate annotation dictionary and return number of genes with more than one CDS
    counter = 0
    for k, v in natsorted(list(indict.items())):
        if len(v["CDS"][0]) > 1:
            counter += 1
    return len(indict), counter


def selectTrainingModels(
    genome, train_dict, tmpdir="/tmp", flank_length=1000, mult_cds_threshold=0.65
):
    """
    Filter and sort gene models from a GFF3 file based on completeness, non-overlapping nature, and number of exons.

    This function processes gene models from a GFF3 file, ensuring they are complete and non-overlapping.
    It sorts the models by the number of exons, preferring those with more exons. The function also
    performs a pairwise DIAMOND search to ensure gene model uniqueness and returns a dictionary of
    filtered and sorted gene models.

    Parameters:
    - genome (str): Path to the genome file in FASTA format.
    - train_dict (dict): A dictionary containing gene models from a GFF3 file.
    - tmpdir (str, optional): Temporary directory for intermediate files (default is "/tmp").
    - flank_length (int, optional): Length of flanking regions to include (default is 1000).
    - mult_cds_threshold (float, optional): Threshold for ratio of multi-CDS genes (default is 0.65).

    Returns:
    - dict: A dictionary of filtered and sorted gene models ready for training.
    """

    def _sortDictCDS(d):
        return len(d[1]["CDS"][0])

    def _sortDict(d):
        return (d[1]["contig"], d[1]["location"][0])

    # setup interlap object
    gene_inter = defaultdict(InterLap)

    # add to InterLap output proteins
    proteins = os.path.join(tmpdir, "training-set.proteins.fa")
    augdmnddb = os.path.join(tmpdir, "training-set.dmnd")
    augblastout = os.path.join(tmpdir, "training-set.self.blast.txt")

    # check number of multi-cds genes
    countGenes, countGenesCDS = count_multi_CDS_genes(train_dict)
    logger.debug(f"{countGenes} training set genes; {countGenesCDS} have multi-CDS")

    # calculate ratio of multi-CDS genes
    multiCDSratio = countGenesCDS / countGenes
    multiCDScheck = False
    if multiCDSratio >= mult_cds_threshold:
        multiCDScheck = True
        logger.debug(
            f"multi-CDS ratio is high ({multiCDSratio:.2f}), filtering out single-CDS genes for training"
        )

    with open(proteins, "w") as protout:
        for k, v in natsorted(list(train_dict.items())):
            # ensure that gene model is complete
            if v["protein"][0].startswith("M") and v["protein"][0].endswith("*"):
                if multiCDScheck:
                    if len(v["CDS"][0]) > 1:
                        # add to interlap object and write protein out
                        gene_inter[v["contig"]].add(
                            (
                                v["location"][0],
                                v["location"][1],
                                v["strand"],
                                k,
                                len(v["CDS"][0]),
                            )
                        )
                        protout.write(
                            f">{k}___{len(v['CDS'][0])}\n{softwrap(v['protein'][0])}\n"
                        )
                else:
                    gene_inter[v["contig"]].add(
                        (
                            v["location"][0],
                            v["location"][1],
                            v["strand"],
                            k,
                            len(v["CDS"][0]),
                        )
                    )
                    protout.write(
                        f">{k}___{len(v['CDS'][0])}\n{softwrap(v['protein'][0])}\n"
                    )

    # make sure gene models are unique, so do pairwise diamond search @ 80% identity
    cmd = ["diamond", "makedb", "--in", proteins, "--db", augdmnddb]
    runSubprocess(cmd, logger, only_failed=True)
    cmd = [
        "diamond",
        "blastp",
        "--query",
        proteins,
        "--db",
        augdmnddb,
        "--more-sensitive",
        "-o",
        augblastout,
        "-f",
        "6",
        "qseqid",
        "sseqid",
        "pident",
        "--query-cover",
        "80",
        "--subject-cover",
        "80",
        "--id",
        "80",
        "--no-self-hits",
    ]
    runSubprocess(cmd, logger, only_failed=True)
    blast_results = []
    with open(augblastout, "r") as blast:
        for line in blast:
            line = line.rstrip()
            line = line.replace("___", "\t")
            blast_results.append(line.split("\t"))
    sortedBlast = natsorted(blast_results, key=lambda x: int(x[1]), reverse=True)
    blastignore = []
    for hit in sortedBlast:
        if hit[0] in blastignore or hit[2] in blastignore:
            continue
        if int(hit[1]) >= int(hit[3]):
            if hit[2] not in blastignore:
                blastignore.append(hit[2])
        else:
            if hit[0] not in blastignore:
                blastignore.append(hit[0])
    logger.debug("{:,} models fail blast identity threshold".format(len(blastignore)))

    GenesPass = {}
    for k, v in natsorted(list(train_dict.items())):
        if k not in blastignore and k not in GenesPass:
            loc = sorted([v["location"][0], v["location"][1]])
            if loc in gene_inter[v["contig"]]:
                hits = list(gene_inter[v["contig"]].find(loc))
                sortedHits = sorted(hits, key=lambda x: int(x[4]), reverse=True)
                validHits = []
                for y in sortedHits:
                    if y[3] not in blastignore and y[3] != k:
                        validHits.append(y)
                if len(validHits) > 0:
                    if validHits[0][3] not in GenesPass:
                        GenesPass[validHits[0][3]] = train_dict.get(validHits[0][3])
                else:
                    GenesPass[k] = v

    # now sort dictionary number of exons if multiCDScheck else just location
    if multiCDScheck:
        sGenes = sorted(iter(GenesPass.items()), key=_sortDictCDS, reverse=True)
    else:
        sGenes = sorted(iter(GenesPass.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    logger.info(
        "{:,} of {:,} models pass training parameters".format(
            len(sortedGenes), len(train_dict)
        )
    )
    # normalize the names, but also write each gene +/- 1 kb to file
    fa = fasta2dict(genome)
    final = {}
    for i, (k, v) in enumerate(natsorted(list(sortedGenes.items()))):
        gene_name = "g_" + str(i + 1)
        v["ids"] = [gene_name + "-T1"]
        final[gene_name] = v
        cstart = min(sorted(v["CDS"][0], key=lambda x: x[0])[0]) - flank_length
        if cstart < 0:
            cstart = 0
        cend = max(sorted(v["CDS"][0], key=lambda x: x[0])[-1]) + flank_length
        if cend > len(fa[v["contig"]]):
            cend = len(fa[v["contig"]])
        gSeq = fa[v["contig"]][cstart:cend]
        assert (cend - cstart) == len(gSeq)
        cdsCoords = [(x[0] - cstart, x[1] - cstart) for x in v["CDS"][0]]
        if v["strand"] == "+":
            cdsStart = min(cdsCoords[0])
            cdsEnd = max(cdsCoords[-1])
        else:
            cdsStart = min(cdsCoords[-1])
            cdsEnd = max(cdsCoords[0])
        assert (len(gSeq) - (cdsEnd - cdsStart)) <= 2000
        final[gene_name]["train_dna"] = gSeq
        final[gene_name]["train_coords"] = cdsCoords

    # return
    return final


def trainmodels2dict(genome, models, flank_length=1000):
    """
    Parse existing gene models and return a modified dictionary with flanking sequences.

    This function processes gene models from a GFF3 file and a genome FASTA file. It extracts
    sequences with specified flanking regions for each gene model, ensuring they are not pseudogenes
    or non-coding RNAs. The function returns a dictionary with the modified gene models, including
    the extracted sequences and coordinates.

    Parameters:
    - genome (str): Path to the input genome file in FASTA format.
    - models (str): Path to the input gene models file in GFF3 format.
    - flank_length (int, optional): Length of flanking sequences to include (default is 1000).

    Returns:
    - dict: A dictionary containing modified gene models with flanking sequences and coordinates.
    """
    # existing models parse and return the modified dictionary
    fa = fasta2dict(genome)
    Genes = gff2dict(models, genome)
    final = {}
    for k, v in Genes.items():
        if v["pseudo"] is True or "ncRNA" in v["type"]:
            continue
        try:
            cstart = min(sorted(v["CDS"][0], key=lambda x: x[0])[0]) - flank_length
        except IndexError:
            print(v)
            raise SystemExit(1)
        if cstart < 0:
            cstart = 0
        cend = max(sorted(v["CDS"][0], key=lambda x: x[0])[-1]) + flank_length
        if cend > len(fa[v["contig"]]):
            cend = len(fa[v["contig"]])
        gSeq = fa[v["contig"]][cstart:cend]
        assert (cend - cstart) == len(gSeq)
        cdsCoords = [(x[0] - cstart, x[1] - cstart) for x in v["CDS"][0]]
        if v["strand"] == "+":
            cdsStart = min(cdsCoords[0])
            cdsEnd = max(cdsCoords[-1])
        else:
            cdsStart = min(cdsCoords[-1])
            cdsEnd = max(cdsCoords[0])
        assert (len(gSeq) - (cdsEnd - cdsStart)) <= 2000
        final[k] = v
        final[k]["train_dna"] = gSeq
        final[k]["train_coords"] = cdsCoords
    return final
