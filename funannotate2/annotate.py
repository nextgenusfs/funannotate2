import json
import os
import shutil
import sys
import time

# from natsort import natsorted
from collections import OrderedDict

from gfftk.convert import _dict2proteins, _dict2transcripts
from gfftk.fasta import fasta2lengths
from gfftk.genbank import dict2tbl, table2asn, tbl2dict
from gfftk.gff import dict2gff3, gff2dict
from gfftk.stats import annotation_stats

from .config import env
from .log import finishLogging, startLogging, system_info
from .name_cleaner import (
    NameCleaner,
    write_new_valid_annotations,
    write_problematic_annotations,
)

# from .fastx import fasta2chunks
from .search import (
    busco2tsv,
    busco_search,
    dbcan2tsv,
    dbcan_search,
    digitize_sequences,
    merops2tsv,
    merops_blast,
    parse_annotations,
    pfam2tsv,
    pfam_search,
    swissprot2tsv,
    swissprot_blast,
)
from .utilities import (
    checkfile,
    choose_best_busco_species,
    create_directories,
    create_tmpdir,
    find_files,
    load_json,
    lookup_taxonomy,
    naming_slug,
    get_odb_version,
    validate_busco_lineage,
)
from .config import busco_taxonomy


def _sortDict(d):
    return (d[1]["location"][0], d[1]["location"][1])


def annotate(args):
    # parse the input and get files that you need
    taxonomy = None
    if args.input_dir:
        if os.path.isdir(args.input_dir):
            if not args.fasta:
                fasta_files = find_files(
                    os.path.join(args.input_dir, "predict_results"), ".fasta"
                )
                if len(fasta_files) == 1:
                    args.fasta = os.path.abspath(fasta_files[0])
            if not args.tbl:
                tbl_files = find_files(
                    os.path.join(args.input_dir, "predict_results"), ".tbl"
                )
                if len(tbl_files) == 1:
                    args.tbl = os.path.abspath(tbl_files[0])
            if not args.gff3:
                gff_files = find_files(
                    os.path.join(args.input_dir, "predict_results"), ".gff3"
                )
                if len(gff_files) == 1:
                    args.gff3 = os.path.abspath(gff_files[0])
            if not args.out:
                args.out = args.input_dir
            if not args.species:
                summaryjson = find_files(
                    os.path.join(args.input_dir, "predict_results"), "summary.json"
                )
                if len(summaryjson) == 1:
                    predict_json = load_json(summaryjson[0])
                    args.species = predict_json["species"]
                    taxonomy = predict_json["taxonomy"]
                    args.strain = predict_json["strain"]
        else:
            sys.stderr.write("ERROR: -i,--input-dir is not a directory\n")
            raise SystemExit(1)
    if not args.fasta:
        sys.stderr.write(
            "ERROR: -f,--fasta is required if not passing an existing funannotate2 --input-dir \n"
        )
        raise SystemExit(1)
    if not args.tbl and not args.gff3:
        sys.stderr.write(
            "ERROR: annotation required [-t,--tbl or -g,--gff3] if not passing an existing funannotate2 --input-dir \n"
        )
        raise SystemExit(1)
    if not args.out:
        sys.stderr.write(
            "ERROR: -o,--out is required if not passing an existing funannotate2 --input-dir \n"
        )
        raise SystemExit(1)

    # create output directories
    misc_dir, res_dir, log_dir = create_directories(args.out, base="annotate")

    # start logger
    logger = startLogging(logfile=os.path.join(log_dir, "funannotate-annotate.log"))
    log = logger.info
    system_info(log)

    # check dependencies, log to logfile

    # log the parsed input files
    if args.input_dir:
        logger.info(
            f"Parsed input files from --input-dir {args.input_dir}\n  --fasta {args.fasta}\n  --tbl {args.tbl}\n  --gff3 {args.gff3}\n  --out {args.out}"
        )

    # create a tmpdir for some files
    tmp_dir = create_tmpdir(args.tmpdir, base="annotate")
    logger.info(f"temporary files located in: {tmp_dir}")

    # we need to get the protein sequences and gene models in Genes object
    if args.tbl:
        Genes, parse_errors = tbl2dict(args.tbl, args.fasta)
        if len(parse_errors) > 1:
            logger.warning("There were errors parsing the TBL annotation file")
    elif args.gff3:
        Genes = gff2dict(args.gff3, args.fasta)

    # Get genome stats from the annotation
    genome_stats_data = annotation_stats(Genes)
    logger.info("Parsed genome stats:")
    logger.info(f"\n{json.dumps(genome_stats_data, indent=2)}")
    Proteins = os.path.join(misc_dir, "proteome.fasta")

    # Write protein sequences
    _dict2proteins(Genes, output=Proteins, strip_stop=True)

    # for pyhmmer load query sequences into digitized list
    # will not use chunked because pyhmmer is plenty fast and parallized natively
    digital_seqs = digitize_sequences(Proteins)

    # run Pfam mapping
    pfam_all = os.path.join(misc_dir, "pfam.results.json")
    pfam_annots = os.path.join(misc_dir, "annotations.pfam.tsv")
    if not checkfile(pfam_annots):
        logger.info("Annotating proteome with pyhmmer against the Pfam-A database")
        start = time.time()
        pfam = pfam_search(digital_seqs, cpus=args.cpus)
        end = time.time()
        logger.info(
            f"Pfam-A search resulted in {len(pfam)} hits and finished in {round(end - start, 2)} seconds"
        )
        pfam_dict = pfam2tsv(pfam, pfam_all, pfam_annots)
    else:
        pfam_dict = parse_annotations(pfam_annots)
        logger.info(
            f"Existing Pfam-A results found, loaded annotations for {len(pfam_dict)} gene models"
        )

    # now lets try dbCAN with same method
    dbcan_all = os.path.join(misc_dir, "dbcan.results.json")
    dbcan_annots = os.path.join(misc_dir, "annotations.dbcan.tsv")
    if not checkfile(dbcan_annots):
        logger.info(
            "Annotating proteome with pyhmmer against the dbCAN (CAZyme) database"
        )
        start = time.time()
        dbcan = dbcan_search(digital_seqs, cpus=args.cpus)
        end = time.time()
        logger.info(
            f"dbCAN search resulted in {len(dbcan)} hits and finished in {round(end - start, 2)} seconds"
        )
        dbcan_dict = dbcan2tsv(dbcan, dbcan_all, dbcan_annots)
    else:
        dbcan_dict = parse_annotations(dbcan_annots)
        logger.info(
            f"Existing dbCAN results found, loaded annotations for {len(dbcan_dict)} gene models"
        )

    # diamond based routines below, also no need to run on the split prots as it parallizes properly
    # now can map to UniProtKB/Swiss-Prot
    swiss_all = os.path.join(misc_dir, "uniprot-swissprot.results.json")
    swiss_annots = os.path.join(misc_dir, "annotations.uniprot-swissprot.tsv")
    if not checkfile(swiss_annots):
        logger.info(
            "Annotating proteome with diamond against the UniProtKB/Swiss-Prot database"
        )
        start = time.time()
        swiss = swissprot_blast(Proteins, cpus=args.cpus)
        end = time.time()
        logger.info(
            f"UniProtKB/Swiss-Prot search resulted in {len(swiss)} hits and finished in {round(end - start, 2)} seconds"
        )
        swiss_dict = swissprot2tsv(swiss, swiss_all, swiss_annots)
    else:
        swiss_dict = parse_annotations(swiss_annots)
        logger.info(
            f"Existing UniProtKB/Swiss-Prot results found, loaded annotations for {len(swiss_dict)} gene models"
        )

    # merops
    merops_all = os.path.join(misc_dir, "merops.results.json")
    merops_annots = os.path.join(misc_dir, "annotations.merops.tsv")
    if not checkfile(merops_annots):
        logger.info(
            "Annotating proteome with diamond against the MEROPS protease database"
        )
        start = time.time()
        merops = merops_blast(Proteins, cpus=args.cpus)
        end = time.time()
        logger.info(
            f"MEROPS search resulted in {len(merops)} hits and finished in {round(end - start, 2)} seconds"
        )
        merops_dict = merops2tsv(merops, merops_all, merops_annots)
    else:
        merops_dict = parse_annotations(merops_annots)
        logger.info(
            f"Existing MEROPS results found, loaded annotations for {len(merops_dict)} gene models"
        )

    # busco proteome analysis
    busco_all = os.path.join(misc_dir, "busco.results.json")
    busco_annots = os.path.join(misc_dir, "annotations.busco.tsv")
    odb_version = get_odb_version(
        os.path.join(os.path.dirname(__file__), "downloads.json")
    )
    if not checkfile(busco_annots):
        if not taxonomy:
            # get taxonomy information
            taxonomy = lookup_taxonomy(args.species)

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
        busco_model_path = os.path.join(
            env["FUNANNOTATE2_DB"], f"{busco_species}_{odb_version}"
        )

        # run busco proteome screen
        logger.info(
            f"BUSCOlite [conserved ortholog] search using {busco_species} models"
        )
        start = time.time()
        busco_results = busco_search(
            Proteins, busco_model_path, cpus=args.cpus, logger=logger
        )
        end = time.time()
        busco_dict = busco2tsv(busco_results, busco_model_path, busco_all, busco_annots)
        logger.info(
            f"BUSCOlite search resulted in {len(busco_dict)} hits and finished in {round(end - start, 2)} seconds"
        )
    else:
        busco_dict = parse_annotations(busco_annots)
        logger.info(
            f"Existing BUSCOlite results found, loaded annotations for {len(busco_dict)} gene models"
        )

    # load any annotations from cli and merge rest of annotations
    all_annotations = [pfam_dict, dbcan_dict, swiss_dict, merops_dict, busco_dict]
    # we need to look for any additional annotations that might be added by f2a
    annots_in_dir = find_files(misc_dir, ".annotations.txt")
    for annotfile in annots_in_dir:
        a = parse_annotations(annotfile)
        logger.info(f"Loaded {len(a)} annotations from {annotfile}")
        all_annotations.append(a)
    if args.annotations:
        for annotfile in args.annotations:
            a = parse_annotations(annotfile)
            logger.info(f"Loaded {len(a)} annotations from {annotfile}")
            all_annotations.append(a)

    # merge annotations into the gene/funannotate dictionary
    merged = {}
    sources = {}
    for annot in all_annotations:
        for g, f in annot.items():
            if g not in merged:
                merged[g] = f
                for key, value in f.items():
                    if key not in sources:
                        sources[key] = 1
                    else:
                        sources[key] += 1
            else:
                for key, value in f.items():
                    if key not in merged[g]:
                        merged[g][key] = value
                    else:
                        merged[g][key] += value
                    if key not in sources:
                        sources[key] = 1
                    else:
                        sources[key] += 1
    logger.info(f"Found functional annotation for {len(merged)} gene models")
    logger.info(f"Annotation sources: {sources}")

    # Clean gene names and product descriptions using the curated database
    logger.info("Cleaning gene names and product descriptions using curated database")
    if args.curated_names:
        logger.info(f"Using custom annotations from {args.curated_names}")
        name_cleaner = NameCleaner(custom_file=args.curated_names)
    else:
        name_cleaner = NameCleaner()

    # Write problematic annotations to file for manual curation
    problematic_file = os.path.join(res_dir, "Gene2Products.need-curating.txt")
    num_problematic = write_problematic_annotations(merged, problematic_file)
    if num_problematic > 0:
        logger.info(
            f"Found {num_problematic} problematic gene names/products that need manual curation"
        )
        logger.info(f"See {problematic_file} for details")

    # Write new valid annotations to file for potential addition to curated database
    new_valid_file = os.path.join(res_dir, "Gene2Products.new-valid.txt")
    num_new_valid = write_new_valid_annotations(merged, new_valid_file)
    if num_new_valid > 0:
        logger.info(
            f"Found {num_new_valid} new valid gene names/products that could be added to the curated database"
        )
        logger.info(f"See {new_valid_file} for details")

    # First, clean the merged annotations
    cleaned_merged = {}
    for gene_id, annot in merged.items():
        # Process the annotation to clean names and products, passing the gene_id for custom annotations
        cleaned_annot = name_cleaner.process_annotation(annot, gene_id=gene_id)
        cleaned_merged[gene_id] = cleaned_annot

    # now loop through the annotation object and add functional annotation
    Annotation = {}
    for k, v in Genes.items():
        n = v.copy()
        for i, x in enumerate(n["ids"]):
            if x in cleaned_merged:  # then functional annotation to add
                fa = cleaned_merged.get(x)
                if "product" in fa:
                    n["product"][i] = fa["product"][0]
                if "db_xref" in fa:
                    n["db_xref"][i] = fa["db_xref"]
                if "name" in fa:
                    n["name"] = fa["name"][0]
                    if len(fa["name"]) > 1:
                        n["gene_synonym"] += fa["name"][1:]
                if "note" in fa:
                    n["note"][i] = fa["note"]
                if "ec_number" in fa:
                    n["ec_number"][i] = fa["ec_number"]
                if "go_terms" in fa:
                    n["go_terms"][i] = fa["go_terms"]
        Annotation[k] = n

    logger.info("Converting to GenBank format")
    # now write TBL outputfile
    finalTBL = os.path.join(res_dir, f"{naming_slug(args.species, args.strain)}.tbl")
    # first sort the dictionary to get order for tbl
    sGenes = sorted(iter(Annotation.items()), key=_sortDict)
    sortedGenes = OrderedDict(sGenes)
    scaff2genes = {}
    for k, v in list(sortedGenes.items()):
        if v["contig"] not in scaff2genes:
            scaff2genes[v["contig"]] = [k]
        else:
            scaff2genes[v["contig"]].append(k)
    # get contig lengths
    scaffLen = fasta2lengths(args.fasta)
    # finally write output
    errors, _, _, _ = dict2tbl(
        sortedGenes,
        scaff2genes,
        scaffLen,
        "CFMR",
        "12345",
        [],
        output=finalTBL,
        annotations=True,
        external=True,
    )
    if len(errors) > 0:
        logger.warning(
            "Errors detected in creation of TBL file:\n{}".format("\n".join(errors))
        )
    # now we want to run table2asn
    finalGBK = os.path.join(res_dir, f"{naming_slug(args.species, args.strain)}.gbk")
    table2asn(
        finalTBL,
        args.fasta,
        organism=args.species,
        strain=args.strain,
        tmpdir=misc_dir,
        table=1,
        cleanup=False,
        output=finalGBK,
    )
    # fetch files from the table2asn directory

    # now write remaining files
    finalGFF3 = os.path.join(res_dir, f"{naming_slug(args.species, args.strain)}.gff3")
    finalProteins = os.path.join(
        res_dir, f"{naming_slug(args.species, args.strain)}.proteins.fa"
    )
    finalTranscripts = os.path.join(
        res_dir, f"{naming_slug(args.species, args.strain)}.transcripts.fa"
    )
    finalFA = os.path.join(res_dir, f"{naming_slug(args.species, args.strain)}.fasta")
    finalSummary = os.path.join(
        res_dir, f"{naming_slug(args.species, args.strain)}.summary.json"
    )

    # Write GFF3 file directly from the annotation object
    logger.info("Writing rest of the output annotation files")
    dict2gff3(sortedGenes, output=finalGFF3)

    # Extract protein sequences directly from the annotation object
    _dict2proteins(sortedGenes, output=finalProteins, strip_stop=True)

    # Extract transcript sequences directly from the annotation object
    _dict2transcripts(sortedGenes, output=finalTranscripts)

    # Copy the input FASTA to the results directory
    shutil.copy2(args.fasta, finalFA)

    # Generate summary statistics directly from the annotation object
    stats = annotation_stats(sortedGenes)
    with open(finalSummary, "w") as f:
        json.dump(stats, f, indent=4)

    # Print summary to log
    logger.info("Annotation Summary:")
    logger.info(f"\n{json.dumps(stats, indent=2)}")

    # finish
    finishLogging(log, vars(sys.modules[__name__])["__name__"])
