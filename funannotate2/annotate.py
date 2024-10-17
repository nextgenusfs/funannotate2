import sys
import os
import json
import shutil
import time
from natsort import natsorted
from .utilities import create_directories, create_tmpdir, find_files
from .log import startLogging, system_info, finishLogging
from .fastx import fasta2chunks
from .search import (
    digitize_sequences,
    pfam_search,
    dbcan_search,
    dbcan2tsv,
    swissprot_blast,
    swissprot2tsv,
    pfam2tsv,
)
from .config import env
from gfftk.gff import gff2dict, dict2gff3
from gfftk.stats import annotation_stats
from gfftk.convert import _dict2proteins, gff2tbl, tbl2gbff
from gfftk.genbank import tbl2dict
from buscolite.busco import runbusco


def annotate(args):
    # parse the input and get files that you need
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

    # we need to get the protein sequences
    if args.tbl:
        Genes, parse_errors = tbl2dict(args.tbl, args.fasta)
        if len(parse_errors) > 1:
            logger.warning("There were errors parsing the TBL annotation file")
    elif args.gff3:
        Genes = gff2dict(args.gff3, args.fasta)

    genome_stats = annotation_stats(Genes)
    logger.info(f"Parsed genome stats:\n{json.dumps(genome_stats, indent=2)}")
    Proteins = os.path.join(misc_dir, "proteome.fasta")
    _dict2proteins(Genes, output=Proteins, strip_stop=True)

    # split input protein files for parallel processing steps
    if args.cpus < 2:
        cpu_count = 2
    else:
        cpu_count = args.cpus
    prot_files = fasta2chunks(
        Proteins,
        cpu_count - 1,
        os.path.join(misc_dir, "prots"),
        prefix="prots_",
        suffix=".fa",
    )
    logger.info(
        f"Gene models were split into {len(prot_files)} chunk(s) for parallel processing."
    )

    # for pyhmmer load query sequences into digitized list
    # will not use chunked because pyhmmer is plenty fast and parallized natively
    digital_seqs = digitize_sequences(Proteins)

    # run Pfam mapping
    logger.info("Annotating proteome with pyhmmer against the Pfam-A database")
    start = time.time()
    pfam = pfam_search(digital_seqs, cpus=args.cpus)
    end = time.time()
    logger.info(
        f"Pfam-A search resulted in {len(pfam)} hits and finished in {round(end-start, 2)} seconds"
    )
    pfam_all = os.path.join(misc_dir, "pfam.results.json")
    pfam_annots = os.path.join(misc_dir, "annotations.pfam.tsv")
    pfam2tsv(pfam, pfam_all, pfam_annots)

    # now lets try dbCAN with same method
    logger.info("Annotating proteome with pyhmmer against the dbCAN (CAZyme) database")
    start = time.time()
    dbcan = dbcan_search(digital_seqs, cpus=args.cpus)
    end = time.time()
    logger.info(
        f"dbCAN search resulted in {len(dbcan)} hits and finished in {round(end-start, 2)} seconds"
    )
    dbcan_all = os.path.join(misc_dir, "dbcan.results.json")
    dbcan_annots = os.path.join(misc_dir, "annotations.dbcan.tsv")
    dbcan2tsv(dbcan, dbcan_all, dbcan_annots)

    # diamond based routines below, also no need to run on the split prots as it parallizes properly
    # now can map to UniProtKB/Swiss-Prot
    logger.info(
        "Annotating proteome with diamond against the UniProtKB/Swiss-Prot database"
    )
    start = time.time()
    swiss = swissprot_blast(Proteins, cpus=args.cpus)
    end = time.time()
    logger.info(
        f"UniProtKB/Swiss-Prot search resulted in {len(swiss)} hits and finished in {round(end-start, 2)} seconds"
    )
    swiss_all = os.path.join(misc_dir, "uniprot-swissprot.results.json")
    swiss_annots = os.path.join(misc_dir, "annotations.uniprot-swissprot.tsv")
    swissprot2tsv(swiss, swiss_all, swiss_annots)

    # merops

    # finish
    finishLogging(log, vars(sys.modules[__name__])["__name__"])
