import sys
import os
import json
import shutil
from collections import OrderedDict
from natsort import natsorted
from .utilities import (
    create_directories,
    create_tmpdir,
    checkfile,
    runProcessJob,
    which_path,
    lookup_taxonomy,
    choose_best_busco_species,
    load_json,
    download,
    runSubprocess,
    naming_slug,
)
from .log import startLogging, system_info, finishLogging
from .fastx import analyzeAssembly, simplify_headers
from .align import align_transcripts, align_proteins
from .evm import evm_consensus
from .abinitio import (
    run_snap,
    run_glimmerhmm,
    run_augustus,
    run_genemark,
    evidence2hints,
    run_trnascan,
)
from .config import env
from gfftk.gff import gff2dict, dict2gff3
from gfftk.stats import annotation_stats
from gfftk.convert import _dict2proteins, gff2tbl, tbl2gbff
from gfftk.consensus import generate_consensus
from buscolite.busco import runbusco


def predict(args):
    # create output directories
    misc_dir, res_dir, log_dir = create_directories(args.out, base="predict")

    # start logger
    logger = startLogging(logfile=os.path.join(log_dir, "funannotate-predict.log"))
    log = logger.info
    system_info(log)

    # check dependencies, log to logfile
    if os.path.isfile(args.params):
        with open(args.params, "r") as infile:
            params = json.load(infile)
    else:
        # here need to load/check if in database
        pass
    logger.info(
        f'Loaded training params for {params["name"]}: {list(params["abinitio"].keys())}'
    )

    # create a tmpdir for some files
    tmp_dir = create_tmpdir(args.tmpdir, base="predict")
    logger.info(f"temporary files located in: {tmp_dir}")

    # load genome and do some QC checks
    logger.info(
        "Loading genome assembly, running QC checks, calculating softmasked regions and assembly gaps"
    )
    GenomeFasta = os.path.join(misc_dir, "genome.fasta")
    contig_name_map = simplify_headers(args.fasta, GenomeFasta)
    maskedRegions = os.path.join(misc_dir, "softmasked-regions.bed")
    asmGaps = os.path.join(misc_dir, "assembly-gaps.bed")
    stats, bad_names, nuc_errors, contigs = analyzeAssembly(
        GenomeFasta,
        maskedRegions,
        asmGaps,
        header_max=args.header_length,
        split=tmp_dir,
        cpus=args.cpus,
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

    # align transcripts if passed
    TranAlign = os.path.join(misc_dir, "transcript-alignments.gff3")
    TranGenes = os.path.join(misc_dir, "predictions.gapmm2-gene.gff3")
    if args.transcripts and not checkfile(TranAlign):
        align_transcripts(
            GenomeFasta,
            args.transcripts,
            TranAlign,
            TranGenes,
            cpus=args.cpus,
            max_intron=args.max_intron,
            log=logger,
        )
    elif checkfile(TranAlign) and checkfile(TranGenes):
        log("Existing transcript alignments found, will re-use and continue")

    # align proteins if passed
    ProtAlign = os.path.join(misc_dir, "protein-alignments.gff3")
    ProtGenes = os.path.join(misc_dir, "predictions.miniprot-gene.gff3")
    if args.proteins and not checkfile(ProtAlign):
        align_proteins(
            GenomeFasta,
            args.proteins,
            ProtAlign,
            ProtGenes,
            cpus=args.cpus,
            max_intron=args.max_intron,
            log=logger,
        )
    elif checkfile(ProtAlign) and checkfile(ProtGenes):
        log("Existing protein alignments found, will re-use and continue")

    # lets see if already run, get list of expected output files
    Consensus = os.path.join(misc_dir, "consensus.predictions.gff3")
    ab_files_missing = []
    for ab in list(params["abinitio"].keys()):
        if not checkfile(os.path.join(misc_dir, f"predictions.{ab}.gff3")):
            ab_files_missing.append(ab)

    if len(ab_files_missing) > 0:  # all or nothing for now
        if "augustus" in params["abinitio"]:  # generate hintsfiles
            logger.info("Parsing alignments and generating hintsfile for augustus")
            evidence2hints(ProtAlign, TranAlign, contigs, tmp_dir)

        # now run ab intio in parallel
        abinit_cmds = []
        for c in contigs:
            abinit_cmds.append((c, params, logger))

        logger.info(f"Running ab initio gene predictions using {args.cpus} cpus")
        runProcessJob(abinitio_wrapper, abinit_cmds, cpus=args.cpus)

        # get all predictions
        gene_counts = {}
        abinitio_preds = []
        if checkfile(ProtGenes):
            abinitio_preds.append(ProtGenes)
        if checkfile(TranGenes):
            abinitio_preds.append(TranGenes)
        for ab in list(params["abinitio"].keys()):
            gene_counts[ab] = 0
            if ab == "augustus":  # split hiq and regular
                gene_counts["augustus-hiq"] = 0
                abinitio_preds.append(os.path.join(misc_dir, f"predictions.{ab}.gff3"))
                abinitio_preds.append(
                    os.path.join(misc_dir, f"predictions.{ab}-hiq.gff3")
                )
                with open(os.path.join(misc_dir, f"predictions.{ab}.gff3"), "w") as aug:
                    with open(
                        os.path.join(misc_dir, f"predictions.{ab}-hiq.gff3"), "w"
                    ) as hiq:
                        aug.write("##gff-version 3\n")
                        hiq.write("##gff-version 3\n")
                        for f in natsorted(os.listdir(tmp_dir)):
                            if f.endswith(f"{ab}.gff3"):
                                with open(os.path.join(tmp_dir, f), "r") as infile:
                                    for line in infile:
                                        if line.startswith("#"):
                                            continue
                                        if "\tgene\t" in line:
                                            if "augustus-hiq" in line:
                                                gene_counts["augustus-hiq"] += 1
                                            else:
                                                gene_counts[ab] += 1
                                        if "augustus-hiq" in line:
                                            hiq.write(line)
                                        else:
                                            aug.write(line)
            else:
                abinitio_preds.append(os.path.join(misc_dir, f"predictions.{ab}.gff3"))
                with open(
                    os.path.join(misc_dir, f"predictions.{ab}.gff3"), "w"
                ) as outfile:
                    outfile.write("##gff-version 3\n")
                    for f in natsorted(os.listdir(tmp_dir)):
                        if f.endswith(f"{ab}.gff3"):
                            with open(os.path.join(tmp_dir, f), "r") as infile:
                                for line in infile:
                                    if line.startswith("#"):
                                        continue
                                    if "\tgene\t" in line:
                                        gene_counts[ab] += 1
                                    outfile.write(line)

        # clean up
        shutil.rmtree(tmp_dir)

        logger.info(
            f"Ab initio predictions finished:\n{json.dumps(gene_counts, indent=2)}"
        )

    else:
        abinitio_preds = []
        for f in os.listdir(misc_dir):
            if f.startswith("predictions") and f.endswith(".gff3"):
                abinitio_preds.append(os.path.join(misc_dir, f))
        logger.info(
            "Existing ab initio predictions found, re-using these files\n - {}".format(
                "\n - ".join(abinitio_preds)
            )
        )
    # ensure we get files in same order
    abinitio_preds = sorted(abinitio_preds)

    # dynamically figure out the weights to use based on training test scores
    # or we could test again here with busco?
    # busco with a higher level taxonomic group
    taxonomy = lookup_taxonomy(args.species)
    if taxonomy is False:  # then offline mode, pull from training data if possible
        logger.warning(
            "Unable to access JGI taxonomy lookup, reverting to taxonomy from training data"
        )
        taxonomy = params["taxonomy"]
    busco_tax = choose_best_busco_species(
        {"superkingdom": taxonomy["superkingdom"], "kingdom": taxonomy["kingdom"]}
    )
    busco_model_path = os.path.join(env["FUNANNOTATE_DB"], f"{busco_tax}_odb10")
    if not os.path.isdir(busco_model_path):
        download_urls = load_json(
            os.path.join(os.path.dirname(__file__), "downloads.json")
        )
        busco_url = download_urls["busco"][busco_tax][0]
        busco_tgz = os.path.join(env["FUNANNOTATE_DB"], os.path.basename(busco_url))
        logger.info(f"Downloading {busco_tax}_odb10 model from {busco_url}")
        download(busco_url, busco_tgz, wget=False)
        if os.path.isfile(busco_tgz):
            runSubprocess(
                ["tar", "-zxf", os.path.basename(busco_tgz)],
                logger,
                cwd=env["FUNANNOTATE_DB"],
            )
            if os.path.isdir(busco_model_path):
                os.remove(busco_tgz)

    # now we can loop through the abinitio predictions and run busco for completion
    # write this to file for re-use if consensus file already present?
    abinitio_scores_file = os.path.join(misc_dir, "abinitio_algorithm_scoring.json")
    if checkfile(abinitio_scores_file):
        logger.info("Existing abinitio scoring file present, will re-use and continue")
        with open(abinitio_scores_file, "r") as infile:
            abinitio_scores = json.load(infile)
    else:
        abinitio_scores = {}
        logger.info(
            "Measuring assembly completeness with buscolite for all ab initio predictions"
        )
        for ap in abinitio_preds:
            ProtPreds = os.path.join(misc_dir, os.path.basename(ap) + ".prots.fa")
            gene_models = gff2dict(ap, GenomeFasta)
            _dict2proteins(gene_models, output=ProtPreds)
            if checkfile(ProtPreds):
                d, m, stats, cfg = runbusco(
                    ProtPreds,
                    busco_model_path,
                    mode="proteins",
                    cpus=args.cpus,
                    logger=logger,
                    verbosity=0,
                )
                cov = (stats["single-copy"] + stats["duplicated"]) / float(
                    stats["total"]
                )
            else:
                cov = 0.00
            # measure completeness for each tool
            ab_initio_tool = os.path.basename(ap).split(".")[1]
            if ab_initio_tool == "augustus-hiq":  # add these just augustus
                ab_initio_tool = "augustus"
            if ab_initio_tool not in abinitio_scores:
                abinitio_scores[ab_initio_tool] = {"busco": cov}
            else:
                abinitio_scores[ab_initio_tool]["busco"] += cov
        # add the params based scoring to the dictionary to calculate weights
        for k, v in params["abinitio"].items():
            if k in abinitio_scores:
                abinitio_scores[k]["train"] = v["train_results"]
        # write it to file
        with open(abinitio_scores_file, "w") as outfile:
            json.dump(abinitio_scores, outfile, indent=2)

    logger.info(
        "ab initio models scoring by algorithm:\n{}".format(
            json.dumps(abinitio_scores, indent=2)
        )
    )

    # now we want to calculate weight estimations based on these results
    weightings = calculate_weights(abinitio_scores, args.weights)
    logger.info(f"Calculated ab initio weights from data: {weightings}")
    consensus_tmp = os.path.join(tmp_dir, "evm")
    if not os.path.isdir(consensus_tmp):
        os.makedirs(consensus_tmp)
    # now run consensus
    if not checkfile(Consensus):
        p_aligns = []
        if checkfile(ProtAlign):
            p_aligns.append(ProtAlign)
        t_aligns = []
        if checkfile(TranAlign):
            t_aligns.append(TranAlign)
        if args.consensus == "evm":
            _ = evm_consensus(
                GenomeFasta,
                abinitio_preds,
                p_aligns,
                t_aligns,
                weightings,
                Consensus,
                log=logger,
                tmpdir=consensus_tmp,
                cpus=args.cpus,
            )
        elif args.consensus == "gfftk":
            _ = generate_consensus(
                GenomeFasta,
                abinitio_preds,
                p_aligns,
                t_aligns,
                weightings,
                Consensus,
                repeats=maskedRegions,
                tiebreakers="weights",
                log=logger.info,
                debug=True,
                max_intron=args.max_intron,
                evidence_derived_models=["miniprot-gene", "gapmm2-gene"],
            )
    else:
        logger.info("Existing consensus predictions found, will re-use and continue")

    # predict tRNA
    trna_predictions = os.path.join(misc_dir, "tRNA.predictions.gff3")
    if not checkfile(trna_predictions):
        logger.info("Predicting tRNA genes")
        run_trnascan(GenomeFasta, trna_predictions, cpus=args.cpus, log=logger)
    else:
        logger.info("Existing tRNA predictions found, will re-use and continue")

    # need to combine consensus predictions and RNA and rename gene models
    finalGFF3 = os.path.join(res_dir, f"{naming_slug(args.species, args.strain)}.gff3")
    finalFA = os.path.join(res_dir, f"{naming_slug(args.species, args.strain)}.fasta")
    finalTBL = os.path.join(res_dir, f"{naming_slug(args.species, args.strain)}.tbl")
    finalGBK = os.path.join(res_dir, f"{naming_slug(args.species, args.strain)}.gbk")
    logger.info(
        f"Merging all gene models, sorting, and renaming using locus_tag={args.name}"
    )
    consensus_models = merge_rename_models(
        [Consensus, trna_predictions],
        GenomeFasta,
        finalGFF3,
        locus_tag=args.name,
        contig_map=contig_name_map,
    )
    # generate renaming output files
    shutil.copyfile(args.fasta, finalFA)
    logger.info("Converting to GenBank format")
    gff2tbl(finalGFF3, args.fasta, output=finalTBL, table=1)
    tbl2gbff(
        finalTBL,
        args.fasta,
        output=finalGBK,
        organism=args.species,
        strain=args.strain,
    )

    # get some stats for user
    consensus_stats = annotation_stats(consensus_models)
    logger.info(
        "Annotation statistics:\n{}".format(json.dumps(consensus_stats, indent=2))
    )
    # we are finished here with coding sequences, lets check completeness
    ConsensusProts = os.path.join(misc_dir, "consensus.proteins.fasta")
    _dict2proteins(consensus_models, output=ConsensusProts)
    logger.info("Measuring assembly completeness with buscolite")
    d, m, stats, cfg = runbusco(
        ConsensusProts,
        busco_model_path,
        mode="proteins",
        cpus=args.cpus,
        logger=logger,
        verbosity=0,
    )
    logger.info(
        "Assembly completeness:\n complete={:} [{:.2%}]\n single-copy={:} [{:.2%}]\n fragmented={:} [{:.2%}]\n duplicated={:} [{:.2%}]\n missing={:} [{:.2%}]\n total={:} [{:.2%}]".format(
            stats["single-copy"] + stats["duplicated"],
            ((stats["single-copy"] + stats["duplicated"]) / float(stats["total"])),
            stats["single-copy"],
            (stats["single-copy"] / float(stats["total"])),
            stats["fragmented"],
            (stats["fragmented"] / float(stats["total"])),
            stats["duplicated"],
            (stats["duplicated"] / float(stats["total"])),
            stats["missing"],
            (stats["missing"] / float(stats["total"])),
            stats["total"],
            (stats["total"] / float(stats["total"])),
        )
    )

    # finish
    finishLogging(log, vars(sys.modules[__name__])["__name__"])


def merge_rename_models(gffList, genome, output, locus_tag="FUN_", contig_map={}):
    def _sortDict(d):
        return (d[1]["contig"], d[1]["location"][0])

    Genes = {}
    for gff in gffList:
        Genes = gff2dict(
            os.path.abspath(gff), os.path.abspath(genome), annotation=Genes
        )

    sGenes = natsorted(iter(Genes.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    renamedGenes = {}
    counter = 1
    locus_tag = locus_tag.rstrip("_")
    for k, v in list(sortedGenes.items()):
        locusTag = locus_tag + "_" + str(counter).zfill(6)
        renamedGenes[locusTag] = v
        newIds = []
        for i in range(0, len(v["ids"])):
            newIds.append("{}-T{}".format(locusTag, i + 1))
        renamedGenes[locusTag]["ids"] = newIds
        # map back contig names if present
        renamedGenes[locusTag]["contig"] = contig_map.get(v["contig"], v["contig"])
        counter += 1
    dict2gff3(renamedGenes, output=output, source="funannotate2")
    return renamedGenes


def abinitio_wrapper(contig, params, logger):
    # run ab initio predictions now
    if "snap" in params["abinitio"]:
        run_snap(
            contig,
            params["abinitio"]["snap"]["location"],
            os.path.join(
                os.path.dirname(contig), f"{os.path.basename(contig)}.snap.gff3"
            ),
            log=logger,
        )
    if "glimmerhmm" in params["abinitio"]:
        run_glimmerhmm(
            contig,
            params["abinitio"]["glimmerhmm"]["location"],
            os.path.join(
                os.path.dirname(contig), f"{os.path.basename(contig)}.glimmerhmm.gff3"
            ),
            log=logger,
        )
    if "genemark" in params["abinitio"]:
        if which_path("gmhmme3"):
            run_genemark(
                contig,
                params["abinitio"]["genemark"]["location"],
                os.path.join(
                    os.path.dirname(contig), f"{os.path.basename(contig)}.genemark.gff3"
                ),
                log=logger,
            )
    if "augustus" in params["abinitio"]:
        hints = False
        hintsfile = os.path.join(
            os.path.dirname(contig),
            f"{os.path.basename(contig)}.hintsfile",
        )
        if os.path.isfile(hintsfile):
            hints = hintsfile
        run_augustus(
            contig,
            params["abinitio"]["augustus"]["species"],
            os.path.join(
                os.path.dirname(contig), f"{os.path.basename(contig)}.augustus.gff3"
            ),
            log=logger,
            hints=hints,
            config_path=params["abinitio"]["augustus"]["location"],
        )


def calculate_weights(scores, cli_weights):
    # calculate scores for consensus prediction
    # user defined weights take precedence
    # generally we will work off of augustus data here
    weights = {}
    # first cli_weights is a list
    usr_weights = {}
    if cli_weights:
        for x in cli_weights:
            if ":" in x:
                tool, weight = x.split(":", 1)
                usr_weights[tool] = int(weight)
    aug_busco = scores["augustus"]["busco"]
    aug_aed = scores["augustus"]["train"]["average_aed"]
    aug_exon = scores["augustus"]["train"]["exon_sensitivity"]
    for k, v in scores.items():
        if k in usr_weights:
            weights[k] = usr_weights.get(k)
        else:
            if (
                "train" not in v
            ):  # these are either user entered or from gapmm2/miniprot
                weights[k] = 1
            else:
                busco_diff = aug_busco - v["busco"]  # higher better
                aed_diff = v["train"]["average_aed"] - aug_aed  # lower better
                exon_diff = aug_exon - v["train"]["exon_sensitivity"]  # higher better
                change = busco_diff + aed_diff + exon_diff
                if change <= -0.20:  # 20% better than augustus
                    weights[k] = 4
                elif change < 0:
                    weights[k] = 3
                elif change == 0:
                    weights[k] = 2
                elif change > 0.20:  # 20% worse than augustus
                    weights[k] = 1
                else:
                    weights[k] = 2
    # add hiq augustus
    weights["augustus-hiq"] = weights["augustus"] + 2
    # add any other usr_weights
    for k, v in usr_weights.items():
        if k not in weights:
            weights[k] = v
    # now we want a list with tool:weight annotations
    weight_list = []
    for k, v in weights.items():
        weight_list.append(f"{k}:{v}")
    return weight_list
