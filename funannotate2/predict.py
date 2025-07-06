import gzip
import json
import os
import shutil
import sys
from collections import OrderedDict

from buscolite.busco import runbusco
from gfftk.consensus import generate_consensus
from gfftk.convert import _dict2proteins, _dict2transcripts, gff2tbl, tbl2gbff
from gfftk.gff import dict2gff3, gff2dict
from gfftk.stats import annotation_stats
from gfftk.utils import check_file_type
from natsort import natsorted

from .abinitio import (
    evidence2hints,
    run_augustus,
    run_genemark,
    run_glimmerhmm,
    run_snap,
    run_trnascan,
)
from .align import align_mito, align_proteins, align_transcripts
from .config import env
from .database import fetch_pretrained_species
from .fastx import (
    analyzeAssembly,
    annotate_fasta,
    mergefasta,
    simplify_headers_drop,
    softmask_fasta,
)
from .log import finishLogging, startLogging, system_info
from .utilities import (
    checkfile,
    choose_best_busco_species,
    create_directories,
    create_tmpdir,
    download,
    find_files,
    load_json,
    lookup_taxonomy,
    naming_slug,
    runProcessJob,
    runSubprocess,
    which_path,
    get_odb_version,
)


def predict(args):
    """
    Predict gene annotations based on the provided input directory and parameters.

    This function parses the inputs, including parameter files and FASTA files, to predict
    gene annotations. It extracts species information from the parameters file and sets the
    output directory accordingly. The function performs various checks and processes, such as
    loading genome data, aligning transcripts and proteins, running ab initio predictions,
    calculating weights, and generating consensus predictions.

    Parameters:
    - args (argparse.Namespace): Parsed command-line arguments.

    Returns:
    - None
    """
    # parse the inputs, need to do this here
    params = None
    if args.input_dir:  # then can parse some stuff
        if not args.params:
            param_files = find_files(
                os.path.join(args.input_dir, "train_results"), ".params.json"
            )
            if len(param_files) == 1:
                args.params = os.path.abspath(param_files[0])
                with open(args.params, "r") as infile:
                    params = json.load(infile)
            else:
                sys.stderr.write(
                    f"WARNING: unable to find params file in {os.path.join(args.input_dir, 'train_results')}\n{os.listdir(os.path.join(args.input_dir, 'train_results'))}"
                )
        if not args.species:  # load from params file
            args.species = params.get("species")
        if not args.out:
            args.out = args.input_dir
        if not args.fasta:
            fasta_files = find_files(
                os.path.join(args.input_dir, "train_results"),
                (".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz"),
            )
            if len(fasta_files) == 1:
                args.fasta = os.path.abspath(fasta_files[0])
    # check if pre-trained is passed, if so then overwrite any auto globbed data
    if args.params:
        if os.path.isfile(args.params):
            with open(args.params, "r") as infile:
                params = json.load(infile)
        else:
            all_pretrained = fetch_pretrained_species()
            if args.params not in all_pretrained:
                sys.stderr.write(
                    "ERROR: --pretrained-species is not in the database. To see valid species: funannotate2 species"
                )
                raise SystemExit(1)
            else:
                params = all_pretrained.get(args.params)

    # now check arguments
    if not args.out:
        sys.stderr.write(
            "ERROR: -o,--out parameter is missing, exiting. To see arguments: funannotate2 predict --help\n"
        )
        raise SystemExit(1)
    if not args.fasta:
        sys.stderr.write(
            "ERROR: -f,--fasta parameter is missing, exiting. To see arguments: funannotate2 predict --help\n"
        )
        raise SystemExit(1)
    if not args.species:
        sys.stderr.write(
            "ERROR: -s,--species parameter is missing, exiting. To see arguments: funannotate2 predict --help\n"
        )
        raise SystemExit(1)
    if not args.params:
        if params is None:
            sys.stderr.write(
                "ERROR: -p,--params parameter is missing, exiting. To see arguments: funannotate2 predict --help\n"
            )
        raise SystemExit(1)

    # create output directories
    misc_dir, res_dir, log_dir = create_directories(args.out, base="predict")

    # start logger
    logger = startLogging(logfile=os.path.join(log_dir, "funannotate-predict.log"))
    log = logger.info
    system_info(log)

    # output parsed options
    if args.input_dir:  # alert user to parsed input data
        logger.info(
            f'Parsed data from --input-dir {args.input_dir}\n  --fasta {args.fasta}\n  --species "{args.species}"\n  --params {args.params}\n  --out {args.out}'
        )

    # check dependencies, log to logfile
    if params is None:
        if os.path.isfile(args.params):
            with open(args.params, "r") as infile:
                params = json.load(infile)
    logger.info(
        f"Loaded training params for {params['name']}: {list(params['abinitio'].keys())}"
    )

    # create a tmpdir for some files
    tmp_dir = create_tmpdir(args.tmpdir, base="predict")
    logger.info(f"temporary files located in: {tmp_dir}")

    # load genome and do some QC checks
    logger.info(
        "Loading genome assembly, running QC checks, searching for mitochondrial contigs, calculating softmasked regions and assembly gaps"
    )
    mito_contigs, _ = align_mito(args.fasta, cpus=args.cpus)
    if mito_contigs is None:
        logger.warning(
            "Mitochondrial refseq database is not installed, unable to filter contigs"
        )
    if mito_contigs:
        logger.info(
            f"Separating {len(mito_contigs)} mitochondrial contig(s) from the nuclear genome, will recombine at the end of predict\n{mito_contigs}"
        )
    GenomeFasta = os.path.join(misc_dir, "genome.fasta")
    GenomeMito = os.path.join(misc_dir, "genome.mito.fasta")
    contig_name_map = simplify_headers_drop(
        args.fasta, GenomeFasta, GenomeMito, drop=list(mito_contigs.keys())
    )
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

    # if no softmasking done, interject and do it with pytantan
    if float(stats["softmasked"].strip("%")) == 0:
        logger.warning(
            "Genome assembly is not softmasked, running pytantan to quickly softmask repeats"
        )
        os.rename(GenomeFasta, GenomeFasta + ".original")
        softmask_fasta(GenomeFasta + ".original", GenomeFasta)
        stats, bad_names, nuc_errors, contigs = analyzeAssembly(
            GenomeFasta,
            maskedRegions,
            asmGaps,
            header_max=args.header_length,
            split=tmp_dir,
            cpus=args.cpus,
        )

    logger.info(f"Genome stats:\n{json.dumps(stats, indent=2)}")

    # align transcripts if passed
    TranAlign = os.path.join(misc_dir, "transcript-alignments.gff3")
    TranGenes = os.path.join(misc_dir, "predictions.gapmm2-gene.gff3")
    Transcripts = os.path.join(misc_dir, "transcripts.to.align.fasta")
    if args.transcripts and not checkfile(Transcripts):
        n_transcripts, tot_transcripts = mergefasta(args.transcripts, Transcripts)
        logger.info(
            f"Parsed {n_transcripts} [out of {tot_transcripts}] transcripts to align for evidence, derived from:\n{'/n'.join(args.transcripts)}"
        )
    if checkfile(Transcripts) and not checkfile(TranAlign):
        align_transcripts(
            GenomeFasta,
            Transcripts,
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
    Proteins = os.path.join(misc_dir, "proteins.to.align.fasta")
    # use funannotate db as default else add to
    uniprot_db = os.path.abspath(
        os.path.join(env.get("FUNANNOTATE2_DB"), "uniprot_sprot.fasta")
    )
    # get full paths to inputs
    if args.proteins:
        args.proteins = [os.path.abspath(x) for x in args.proteins]
    else:
        args.proteins = []
    if checkfile(uniprot_db) and uniprot_db not in args.proteins:
        args.proteins.append(uniprot_db)
    # merge and dereplicate
    n_prots, tot_prots = mergefasta(args.proteins, Proteins)
    logger.info(
        f"Parsed {n_prots} [out of {tot_prots}] proteins to align for evidence, derived from:\n{'/n'.join(args.proteins)}"
    )
    if checkfile(Proteins) and not checkfile(ProtAlign):
        align_proteins(
            GenomeFasta,
            Proteins,
            ProtAlign,
            ProtGenes,
            cpus=args.cpus,
            max_intron=args.max_intron,
            log=logger,
        )
    elif checkfile(ProtAlign) and checkfile(ProtGenes):
        log("Existing protein alignments found, will re-use and continue")

    # if external predictions passed, process them here
    external_abinitio = []
    external_counts = {}
    if args.external:
        for ext in args.external:
            if checkfile(ext):
                models, sources = sanitize_external(ext, args.fasta, contig_name_map)
                source = list(sources)[0].lower()
                if source in params["abinitio"].keys():
                    source = f"{source}-ext"
                extPred = os.path.join(misc_dir, f"predictions.{source}.gff3")
                dict2gff3(models, output=extPred, debug=False)
                logger.info(
                    f"Parsed {len(models)} external gene models [source={source}] from {ext}"
                )
                external_abinitio.append(extPred)
                external_counts[source] = len(models)
            else:
                logger.error(f"ERROR: -e {ext} is not a valid file, skipping.")

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
        gene_counts = external_counts
        abinitio_preds = external_abinitio
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
    # pull the latest odb version from downloads link
    odb_version = get_odb_version(
        os.path.join(os.path.dirname(__file__), "downloads.json")
    )
    busco_model_path = os.path.join(
        env["FUNANNOTATE2_DB"], f"{busco_tax}_{odb_version}"
    )
    if not os.path.isdir(busco_model_path):
        download_urls = load_json(
            os.path.join(os.path.dirname(__file__), "downloads.json")
        )
        busco_url = download_urls["busco"][busco_tax][0]
        busco_tgz = os.path.join(env["FUNANNOTATE2_DB"], os.path.basename(busco_url))
        logger.info(f"Downloading {busco_tax}_{odb_version} model from {busco_url}")
        download(busco_url, busco_tgz, wget=False)
        if os.path.isfile(busco_tgz):
            runSubprocess(
                ["tar", "-zxf", os.path.basename(busco_tgz)],
                logger,
                cwd=env["FUNANNOTATE2_DB"],
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
            f"Measuring assembly completeness with buscolite [lineage={os.path.basename(busco_model_path)}] for all ab initio predictions"
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
        # run gfftk consensus
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
            debug=False,
            max_intron=args.max_intron,
            evidence_derived_models=["miniprot-gene", "gapmm2-gene"],
            num_processes=args.cpus,
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
    finalSummary = os.path.join(
        res_dir, f"{naming_slug(args.species, args.strain)}.summary.json"
    )
    finalProteins = os.path.join(
        res_dir, f"{naming_slug(args.species, args.strain)}.proteins.fa"
    )
    finalTranscripts = os.path.join(
        res_dir, f"{naming_slug(args.species, args.strain)}.transcripts.fa"
    )
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
    # if there are mitoContigs, add them back but annotate the contig header
    if mito_contigs:
        annotate_fasta(args.fasta, finalFA, ids=list(mito_contigs.keys()))
    else:
        if check_file_type(args.fasta) == "gzipped binary":
            with gzip.open(args.fasta, "rt") as infile, open(finalFA, "w") as outfile:
                for line in infile:
                    outfile.write(line)
        else:
            shutil.copyfile(args.fasta, finalFA)
    logger.info("Converting to GenBank format")
    gff2tbl(finalGFF3, finalFA, output=finalTBL, table=1)
    tbl2gbff(
        finalTBL,
        finalFA,
        output=finalGBK,
        organism=args.species,
        strain=args.strain,
        table=1,
    )
    # write remaining output files
    _dict2proteins(consensus_models, output=finalProteins, strip_stop=True)
    _dict2transcripts(consensus_models, output=finalTranscripts)

    # get some stats for user
    consensus_stats = annotation_stats(consensus_models)
    logger.info(
        "Annotation statistics:\n{}".format(json.dumps(consensus_stats, indent=2))
    )
    # we are finished here with coding sequences, lets check completeness
    logger.info(
        f"Measuring assembly completeness with buscolite [lineage={os.path.basename(busco_model_path)}]"
    )
    d, m, stats, cfg = runbusco(
        finalProteins,
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

    # write a summary json file that can be used downstream
    summary = {
        "genome": args.fasta,
        "taxonomy": taxonomy,
        "species": args.species,
        "strain": args.strain,
        "prediction-params": args.params if args.params is not None else None,
        "busco-models": busco_model_path,
        "busco-completeness": stats,
        "stats": consensus_stats,
    }
    with open(finalSummary, "w") as jsonout:
        json.dump(summary, jsonout, indent=2)

    # finish
    finishLogging(log, vars(sys.modules[__name__])["__name__"])


def sanitize_external(gff3, fasta, contig_map, debug=False):
    """
    Sanitize external gene annotations by swapping contig names based on a provided map.

    This function loads gene annotations from a GFF3 file into a dictionary using the `gff2dict` function. It then swaps contig names using a provided mapping and returns the updated annotations. The function also collects unique gene sources during the process.

    Parameters:
        gff3 (str): Path to the input GFF3 file.
        fasta (str): Path to the corresponding FASTA file.
        contig_map (dict): Mapping of original contig names to new contig names.
        debug (bool, optional): Enable debug mode. Defaults to False.

    Returns:
        tuple: A tuple containing the sanitized gene annotations dictionary and a set of unique sources.
    """
    # load into dicitonary
    Genes = gff2dict(gff3, fasta, debug=debug)
    # swap names of contigs based on the map
    inv_map = {v: k for k, v in contig_map.items()}
    Swap = {}
    source = set()
    for k, v in Genes.items():
        v["contig"] = inv_map.get(v["contig"])
        v["source"] = v["source"].lower()
        Swap[k] = v
        source.add(v["source"])
    return Swap, source


def merge_rename_models(gffList, genome, output, locus_tag="FUN_", contig_map={}):
    """
    Merge and rename gene models from multiple GFF files.

    This function processes a list of GFF files, merging gene models and renaming them
    with a specified locus tag. The gene models are sorted by contig and location, and
    the renamed models are saved in GFF3 format. Optionally, contig names can be mapped
    using a provided dictionary.

    Parameters:
    - gffList (list): List of GFF files to merge gene models from.
    - genome (str): Path to the genome file.
    - output (str): Path to save the merged and renamed gene models in GFF3 format.
    - locus_tag (str, optional): Prefix for the locus tag (default is "FUN_").
    - contig_map (dict, optional): Dictionary mapping contig names (default is {}).

    Returns:
    - dict: Dictionary containing the merged and renamed gene models.
    """

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
    """
    Run ab initio predictions for a given contig based on the specified parameters.

    This function executes various ab initio gene prediction tools (SNAP, GlimmerHMM, GeneMark, and Augustus)
    on the provided contig file. It checks for the presence of each tool in the parameters and runs the
    corresponding prediction, saving the output in GFF3 format. Augustus predictions can optionally use
    hints if a hints file is available.

    Parameters:
    - contig (str): Path to the contig file.
    - params (dict): Dictionary containing ab initio prediction parameters.
    - logger: Logger object for logging messages.

    Returns:
    - None
    """
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
    """
    Calculate weights for consensus prediction based on input scores and user-defined weights.

    This function computes weights for various prediction tools by comparing their scores
    against a reference (typically Augustus). User-defined weights take precedence over
    calculated weights. The function evaluates differences in BUSCO scores, average AED,
    and exon sensitivity to determine the weight for each tool. It also adjusts weights
    for "augustus-hiq" and incorporates any additional user-defined weights.

    Parameters:
    - scores (dict): A dictionary containing scores for different prediction tools.
    - cli_weights (list): A list of user-defined weights in the format "tool:weight".

    Returns:
    - list: A list of strings representing tool-weight pairs in the format "tool:weight".
    """
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
                # can only do augustus busco comparison here
                if v["busco"] >= 0.99:
                    weights[k] = 2
                else:
                    busco_diff = (aug_busco - v["busco"]) * 100
                    # these are percents now
                    if busco_diff <= 1.5:
                        weights[k] = 2
                    else:
                        weights[k] = 1
            else:
                # short circut this logic if busco results are good
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
                    # one more clarifiyer here, interject if busco score < 0.98
                    if v["busco"] >= 0.98:
                        weights[k] = 2
                    else:
                        weights[k] = 1
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
