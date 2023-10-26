import sys
import os
import json
import datetime
import getpass
import random
from collections import defaultdict, OrderedDict
from natsort import natsorted
from buscolite.busco import runbusco
from buscolite.gff import gffwriter
from buscolite.utilities import summary_writer
from gfftk.gff import gff2dict, dict2gff3
from gfftk.fasta import softwrap, fasta2dict
from .interlap import InterLap
from .log import startLogging, system_info, finishLogging
from .fastx import analyzeAssemblySimple
from .config import env
from .utilities import (
    lookup_taxonomy,
    choose_best_augustus_species,
    choose_best_busco_species,
    create_directories,
    checkfile,
    download,
    load_json,
    runSubprocess,
    naming_slug,
    which_path,
)
from .abinitio import train_snap, train_augustus, train_glimmerhmm, train_genemark


def train(args):
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

    # get taxonomy information
    taxonomy = lookup_taxonomy(args.species)
    logger.info(f"Getting taxonomy information\n{json.dumps(taxonomy, indent=2)}")

    # choose best augustus species based on taxonomy
    aug_species = choose_best_augustus_species(taxonomy)
    logger.info(f"Choosing best augustus species based on taxonomy: {aug_species}")

    # choose best busco species
    busco_species = choose_best_busco_species(taxonomy)
    busco_model_path = os.path.join(env["FUNANNOTATE_DB"], f"{busco_species}_odb10")

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
                        env["FUNANNOTATE_DB"], os.path.basename(busco_url)
                    )
                    logger.info(
                        f"Downloading {busco_species}_odb10 model from {busco_url}"
                    )
                    download(busco_url, busco_tgz, wget=False)
                    if os.path.isfile(busco_tgz):
                        runSubprocess(
                            ["tar", "-zxf", os.path.basename(busco_tgz)],
                            logger,
                            cwd=env["FUNANNOTATE_DB"],
                        )
                        if os.path.isdir(busco_model_path):
                            os.remove(busco_tgz)
                log("Running buscolite to generate training set")
                buscolite(
                    args.fasta,
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
        # load  GFF3 training set, load with gfftk
        train_set = gff2dict(args.training_set, args.fasta)
        logger.info(
            f"Training set [{args.training_set}] loaded with {len(train_set)} gene models"
        )
        if len(train_set) == 0:
            logger.critical(
                f"No gene models found in training set: {args.training_set}"
            )
            raise SystemExit(1)

        # now filter training set
        models4training = selectTrainingModels(args.fasta, train_set, tmpdir=misc_dir)
        dict2gff3(models4training, filt_train_models)
    else:
        logger.info(f"Using existing training set: {filt_train_models}")
        models4training = trainmodels2dict(args.fasta, filt_train_models)

    # split into test/train sets
    n_test = int(len(models4training) * 0.20)
    if n_test > 200:
        n_test = 200
    logger.info(
        f"{len(models4training)} gene models selected for training, now splitting into test [n={n_test}] and train [n={len(models4training)-n_test}]"
    )

    test_model_keys = random.sample(list(models4training.keys()), n_test)
    train_model_keys = [x for x in models4training.keys() if x not in test_model_keys]
    test_models = {k: models4training[k] for k in test_model_keys}
    train_models = {k: models4training[k] for k in train_model_keys}
    filt_train_models_final = os.path.join(misc_dir, "training-models.train.gff3")
    dict2gff3(train_models, filt_train_models_final)

    # run augustus training functions
    logger.info(f"Training augustus using training set")
    augustus_train = train_augustus(
        args.fasta,
        train_models,
        test_models,
        folder=misc_dir,
        log=logger,
        cpus=args.cpus,
        optimize=False,
    )
    augustus_train["training_set"] = filt_train_models_final

    # run snap training functions
    logger.info(f"Training snap using training set")
    snap_train = train_snap(
        args.fasta, train_models, test_models, folder=misc_dir, log=logger
    )
    snap_train["training_set"] = filt_train_models_final

    # run glimmerHMM training functions
    logger.info(f"Training glimmerHMM using training set")
    glimm_train = train_glimmerhmm(
        args.fasta,
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
        logger.info(f"Training GeneMark-ES using self-training")
        fungus_flag = False
        if taxonomy["kingdom"] == "Fungi":
            fungus_flag = True
        genemark_train = train_genemark(
            args.fasta,
            train_models,
            test_models,
            folder=misc_dir,
            fungus=fungus_flag,
            cpus=args.cpus,
            log=logger,
        )
        genemark_train["training_set"] = "self training"
        train_data["genemark"] = genemark_train

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
        f"The params.json file can be passed to funannotate2 predict or installed globally with funannotate2 install"
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


def getTrainResults(input):
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


def count_multi_CDS_genes(input):
    # take funannotate annotation dictionary and return number of genes with more than one CDS
    counter = 0
    keepers = []
    for k, v in natsorted(list(input.items())):
        if len(v["CDS"][0]) > 1:
            counter += 1
    return len(input), counter


def selectTrainingModels(genome, train_dict, tmpdir="/tmp", flank_length=1000):
    """
    function to take a GFF3 file and filter the gene models so they complete, non-overalpping,
    and also sort the models by number of exons, the more the better.
    """

    def _sortDict(d):
        return len(d[1]["CDS"][0])

    # setup interlap object
    gene_inter = defaultdict(InterLap)

    # add to InterLap output proteins
    proteins = os.path.join(tmpdir, "training-set.proteins.fa")
    augdmnddb = os.path.join(tmpdir, "training-set.dmnd")
    augblastout = os.path.join(tmpdir, "training-set.self.blast.txt")

    # check number of multi-cds genes
    countGenes, countGenesCDS = count_multi_CDS_genes(train_dict)
    logger.debug(f"{countGenes} training set genes; {countGenesCDS} have multi-CDS")

    multiCDScheck = False
    if countGenesCDS >= 20000:
        multiCDScheck = True

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
                            f'>{k}___{len(v["CDS"][0])}\n{softwrap(v["protein"][0])}\n'
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
                        f'>{k}___{len(v["CDS"][0])}\n{softwrap(v["protein"][0])}\n'
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
            if not hit[2] in blastignore:
                blastignore.append(hit[2])
        else:
            if not hit[0] in blastignore:
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
                    if not y[3] in blastignore and y[3] != k:
                        validHits.append(y)
                if len(validHits) > 0:
                    if not validHits[0][3] in GenesPass:
                        GenesPass[validHits[0][3]] = train_dict.get(validHits[0][3])
                else:
                    GenesPass[k] = v

    # now sort dictionary number of exons
    sGenes = sorted(iter(GenesPass.items()), key=_sortDict, reverse=True)
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
