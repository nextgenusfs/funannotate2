import os
import shutil
import sys
import uuid
from tempfile import NamedTemporaryFile as NTF

import mappy as mp

from .fastx import softwrap
from .log import finishLogging, startLogging, system_info
from .utilities import runThreadJob


def is_duplicated(data, index, tmpdir, min_pident, min_cov):
    """
    Check if the sequence at the given index is duplicated in larger contigs.

    Args:
        data (list): List of dictionaries containing sequence information.
        index (int): Index of the sequence to check for duplication.
        tmpdir (str): Temporary directory path for storing intermediate files.
        min_pident (float): Minimum percentage identity for a match to be considered duplicated.
        min_cov (float): Minimum coverage percentage for a match to be considered duplicated.

    Returns:
        tuple: A tuple containing the header of the sequence, duplication status, percentage identity,
               coverage percentage, and length of the sequence.
    """
    result = False
    pident = 0
    coverage = 0
    # given the genome datatype and index, see if index sequence duplicated in all larger contigs
    with NTF(mode="w", suffix=".fa", prefix="ref_", dir=tmpdir, delete=False) as ref_fa:
        reftmp = ref_fa.name
        for r in data[index + 1 :]:
            ref_fa.write(">{}\n{}\n".format(r["header"], r["sequence"]))
    # generate mappy aligner object
    a = mp.Aligner(reftmp, preset="asm5")
    for hit in a.map(data[index]["sequence"]):
        pident = float(hit.mlen) / int(hit.blen) * 100
        coverage = float(hit.blen) / int(hit.ctg_len) * 100
        if pident > min_pident and coverage > min_cov:
            result = True
            break
    # clean up the tmpfile
    os.remove(reftmp)
    return (data[index]["header"], result, pident, coverage, data[index]["length"])


def load_genome(fafile):
    """
    Parse the genome from the specified FASTA file and calculate the N50 value.

    Args:
        fafile (str): Path to the input FASTA file.

    Returns:
        tuple: A tuple containing a list of dictionaries with headers, lengths, sequences, and statuses sorted by length,
               and the N50 value calculated based on the genome lengths.
    """
    # parse genome, sort by shortest to longest length
    data = []
    for r in mp.fastx_read(fafile, read_comment=False):
        data.append({"header": r[0], "length": len(r[1]), "sequence": r[1], "status": None})
    sortdata = sorted(data, key=lambda x: x["length"])
    lengths = [x["length"] for x in sortdata]

    # Handle empty or very small genomes
    if not lengths:
        return sortdata, 0

    # Calculate N50
    nlist = []
    for x in lengths:
        nlist += [x] * x

    if not nlist:
        return sortdata, 0

    if len(nlist) % 2 == 0 and len(nlist) >= 2:
        medianpos = int(len(nlist) / 2)
        N50 = int((nlist[medianpos] + nlist[medianpos - 1]) / 2)
    elif len(nlist) >= 1:
        medianpos = int(len(nlist) / 2)
        N50 = int(nlist[medianpos])
    else:
        N50 = 0

    return sortdata, N50


def clean(args):
    """
    Process the input arguments, load the genome from a FASTA file, filter contigs based on length,
    check for duplication, rename contigs if specified, and write the filtered contigs to an output file.

    Args:
        args (argparse.Namespace): The parsed command-line arguments.

    Returns:
        None
    """
    logger = startLogging(logfile=args.logfile)
    log = logger.info
    system_info(log)

    # write in a tmpdir
    if args.tmpdir:
        tmpdir = args.tmpdir
    else:
        tmpdir = "clean_{}".format(str(uuid.uuid4()))
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)

    # load genome and sort by length
    genome, n50 = load_genome(args.fasta)

    # filter to minimum length
    genome_ms = [x for x in genome if x["length"] > args.minlen]
    log(
        "Loaded {} contigs; {} are larger than {}; N50 is {} bp".format(
            len(genome), len(genome_ms), args.minlen, n50
        )
    )

    # Check if we have any contigs that meet the minimum length requirement
    if not genome_ms:
        log(
            "No contigs found that meet the minimum length requirement of {} bp".format(args.minlen)
        )
        # Write an empty output file
        with open(args.out, "w") as outfile:
            pass
        log("Wrote 0 contigs to {}".format(args.out))
        shutil.rmtree(tmpdir)
        finishLogging(log, vars(sys.modules[__name__])["__name__"])
        return

    max_idx = len(genome_ms)
    if not args.exhaustive:
        for i, x in enumerate(genome_ms):
            if x["length"] > n50:
                max_idx = i - 1
                break

    # Make sure max_idx is at least 1
    max_idx = max(1, max_idx)

    # now loop through data with thread pool
    log(
        "Checking {} contigs for duplication [minlen={} maxlen={}]".format(
            max_idx, genome_ms[0]["length"], genome_ms[max_idx - 1]["length"]
        )
    )

    # get a list of lists with arguments
    job_arguments = []
    if max_idx > 0:
        for x in range(0, max_idx):
            job_arguments.append([genome_ms, x, tmpdir, args.pident, args.cov])

    # run this will threadpool if we have jobs to run
    results = []
    if job_arguments:
        results = runThreadJob(is_duplicated, job_arguments, cpus=args.cpus, progress=False)

    # parse results
    duplicated = []
    for r in results:
        contig, duplic, pident, cov, clen = r.result()
        if duplic:
            log(
                "{} is duplicated; pident={} coverage={} length={}".format(
                    contig, pident, cov, clen
                )
            )
            duplicated.append(contig)
    if len(duplicated) == 0:
        log("0 Duplicated contigs found")

    # message user if rename passed
    if args.rename:
        log(
            'Renaming contigs using "{}" basename and auto-increment from largest to smallest'.format(
                args.rename
            )
        )

    # write output file
    good_count = 0
    with open(args.out, "w") as outfile:
        for r in reversed(genome_ms):
            if r["header"] not in duplicated:
                good_count += 1
                if args.rename:
                    title = "{}{}".format(args.rename, good_count)
                else:
                    title = r["header"]
                outfile.write(">{}\n{}\n".format(title, softwrap(r["sequence"])))
    log("Wrote {} contigs to {}".format(good_count, args.out))
    shutil.rmtree(tmpdir)
    finishLogging(log, vars(sys.modules[__name__])["__name__"])
