import math
import multiprocessing
import os

import pyfastx
from natsort import natsorted
from pytantan.lib import RepeatFinder, default_scoring_matrix


def softmask_fasta(
    fasta_file,
    softmasked,
    repeat_start=0.005,
    repeat_end=0.05,
    decay=0.9,
    mask=None,
    threshold=0.5,
    protein=False,
    match_score=None,
    mismatch_cost=None,
):
    """
    Softmask a FASTA file by masking repeats based on specified parameters.

    This function reads a FASTA file, identifies repeat sequences using the specified parameters, and writes the softmasked sequences to an output file. The repeats are masked using a specified character or left unchanged if no mask is provided.

    Args:
        fasta_file (str): The input FASTA file to be softmasked.
        softmasked (str): The output file to write the softmasked sequences.
        repeat_start (float, optional): The start threshold for identifying repeats. Defaults to 0.005.
        repeat_end (float, optional): The end threshold for identifying repeats. Defaults to 0.05.
        decay (float, optional): The decay factor for repeat identification. Defaults to 0.9.
        mask (str, optional): The mask character to use. Defaults to None.
        threshold (float, optional): The threshold for masking repeats. Defaults to 0.5.
        protein (bool, optional): Flag indicating if protein sequences are used. Defaults to False.
        match_score (int, optional): The score for a match in the scoring matrix. Defaults to None.
        mismatch_cost (int, optional): The cost for a mismatch in the scoring matrix. Defaults to None.

    Returns:
        None
    """
    # set default scoring matrix
    matrix = default_scoring_matrix(protein, match_score, mismatch_cost)
    repeat_finder = RepeatFinder(
        matrix,
        repeat_start=repeat_start,
        repeat_end=repeat_end,
        decay=decay,
        protein=protein,
    )
    with open(softmasked, "w") as outfile:
        for title, seq in pyfastx.Fasta(fasta_file, build_index=False, full_name=True):
            masked = repeat_finder.mask_repeats(seq, threshold=threshold, mask=mask)
            outfile.write(f">{title}\n{softwrap(masked)}\n")


def countfasta(fasta_file):
    """
    Count the number of sequences in a FASTA file.

    This function reads a FASTA file using the `pyfastx` library and counts the total number of sequences present.

    Args:
        fasta_file (str): Path to the input FASTA file.

    Returns:
        int: Total number of sequences in the FASTA file.
    """
    count = 0
    for seq in pyfastx.Fasta(fasta_file, build_index=False):
        count += 1
    return count


def softwrap(string, every=80):
    """
    Wrap a string to a specified length by inserting newlines at every given interval.

    This function processes the input string and inserts newlines at specified intervals to ensure that each line does not exceed the given length.

    Args:
        string (str): The input string to be wrapped.
        every (int, optional): The maximum length of each line before inserting a newline. Defaults to 80.

    Returns:
        str: The wrapped string with newlines inserted at the specified interval.
    """
    lines = []
    for i in range(0, len(string), every):
        lines.append(string[i : i + every])
    return "\n".join(lines)


def fasta2dict(fasta_file):
    """
    Extract sequences from a FASTA file and store them in a dictionary.

    This function reads a FASTA file using the `pyfastx` library and creates a dictionary where each key is a sequence title and each value is the corresponding sequence.

    Args:
        fasta_file (str): The path to the FASTA file.

    Returns:
        dict: A dictionary where the keys are sequence titles and the values are sequences.
    """
    fa = {}
    for title, seq in pyfastx.Fasta(fasta_file, build_index=False):
        fa[title] = seq
    return fa


def mergefasta(fasta_files, outfile):
    """
    Merge and dereplicate a list of FASTA files.

    This function processes multiple FASTA files, removes duplicate sequences, and writes the unique sequences to an output file. It keeps track of the total number of sequences and the number of unique sequences.

    Args:
        fasta_files (list): List of input FASTA files to merge.
        outfile (str): Output file to write the merged FASTA sequences.

    Returns:
        tuple: A tuple containing the count of unique sequences and total sequences in the merged output.
    """
    # take a list of fasta files, dereplicate, merge
    seen = set()
    n = 0
    o = 0
    with open(outfile, "w") as output:
        for fa in fasta_files:
            for title, seq in pyfastx.Fasta(fa, build_index=False):
                n += 1
                if seq not in seen:
                    output.write(f">{title}\n{seq}\n")
                    seen.add(seq)
                    o += 1
    return o, n


def annotate_fasta(
    fasta_file, outfile, ids=[], annotation="[mcode=4] [location=mitochondrion]"
):
    """
    Annotate specific sequences in a FASTA file with custom information.

    This function opens a FASTA file, annotates sequences with custom information based on provided IDs,
    and writes the annotated sequences to an output file. It uses the `pyfastx` library for reading FASTA files
    and the `softwrap` function to format the sequences.

    Args:
        fasta_file (str): Path to the input FASTA file.
        outfile (str): Path to the output file where annotated sequences will be written.
        ids (list, optional): List of sequence IDs to be annotated. Defaults to an empty list.
        annotation (str, optional): Custom annotation to be added to the specified sequences. Defaults to "[mcode=4] [location=mitochondrion]".

    Returns:
        None
    """
    with open(outfile, "w") as output:
        for title, seq in pyfastx.Fasta(fasta_file, build_index=False):
            if title in ids:
                title = f"{title} {annotation}"
            output.write(f">{title}\n{softwrap(seq)}\n")


def fasta2chunks(fasta_file, chunks, outputdir, prefix="prots_", suffix=".fa"):
    """
    Split a FASTA file into chunks and save each chunk as a separate file in the specified output directory.

    This function divides the input FASTA file into a specified number of chunks, writes each chunk to a separate file,
    and saves these files in the given output directory. It ensures that the output directory exists and uses the
    `pyfastx` library for efficient FASTA file handling.

    Args:
        fasta_file (str): Path to the input FASTA file.
        chunks (int): Number of chunks to split the FASTA file into.
        outputdir (str): Path to the output directory where the chunked files will be saved.
        prefix (str, optional): Prefix to be added to the names of the output files. Defaults to "prots_".
        suffix (str, optional): Suffix to be added to the names of the output files. Defaults to ".fa".

    Returns:
        list: A list of strings, each representing the path to a saved chunked file.
    """
    # list to store the output files
    files = []
    # check if output exists
    if not os.path.isdir(outputdir):
        os.makedirs(outputdir)
    fa = pyfastx.Fasta(fasta_file)
    n_per_chunk = math.ceil(len(fa) / chunks)
    # built the idx ranges to split
    steps = []
    for i, x in enumerate(range(0, len(fa), n_per_chunk)):
        steps.append((i, x, x + n_per_chunk))
    for step in steps:
        outname = os.path.join(outputdir, f"{prefix}{step[0] + 1}{suffix}")
        files.append(outname)
        with open(outname, "w") as outfile:
            for idx in range(step[1], step[2]):
                try:
                    outfile.write(fa[idx].raw)
                except IndexError:
                    continue
    return files


def simplify_headers(inputfile, outputfile, base="contig_"):
    """
    Simplify the headers of contigs in the input FASTA file and write the simplified sequences to an output file.

    This function reads the input FASTA file, simplifies the headers by assigning new names based on the specified base and order, and writes the simplified sequences to the output file. It maintains the original order of contigs.

    Args:
        inputfile (str): Path to the input FASTA file.
        outputfile (str): Path to the output file where the simplified sequences will be written.
        base (str, optional): Base name for the simplified contigs. Defaults to "contig_".

    Returns:
        dict: A dictionary mapping the new simplified contig names to the original headers.
    """
    # will keep the same order of contigs, just simplify the headers
    names = {}
    with open(outputfile, "w") as outfile:
        for i, (title, seq) in enumerate(pyfastx.Fasta(inputfile, build_index=False)):
            names[f"{base}{i + 1}"] = title
            outfile.write(f">{base}{i + 1}\n{softwrap(seq)}\n")
    return names


def simplify_headers_drop(inputfile, keepfile, dropfile, base="contig_", drop=[]):
    """
    Simplify headers of contigs in the input FASTA file while maintaining order.

    This function simplifies the headers of contigs in the input FASTA file by replacing them with a base string followed by a numerical index.
    Headers specified in the 'drop' list are written to the 'dropfile' with wrapped sequences, while the simplified headers and sequences are written to the 'keepfile'.

    Args:
        inputfile (str): Path to the input FASTA file.
        keepfile (str): Path to the output file to store simplified headers and sequences of contigs to be kept.
        dropfile (str): Path to the output file to store headers and sequences of contigs to be dropped.
        base (str, optional): Base string to be used for simplifying headers. Defaults to "contig_".
        drop (list, optional): List of headers to be dropped. Defaults to an empty list.

    Returns:
        dict: A dictionary mapping simplified headers to original headers.
    """
    # will keep the same order of contigs, just simplify the headers
    names = {}
    with open(keepfile, "w") as outfile:
        with open(dropfile, "w") as dropout:
            for i, (title, seq) in enumerate(
                pyfastx.Fasta(inputfile, build_index=False)
            ):
                if title in drop:
                    dropout.write(f">{title}\n{softwrap(seq)}\n")
                else:
                    names[f"{base}{i + 1}"] = title
                    outfile.write(f">{base}{i + 1}\n{softwrap(seq)}\n")
    return names


def list2groups(L):
    """
    Identify groups of continuous numbers in a list.

    This function processes a list of numbers and identifies groups of consecutive numbers. It yields a tuple representing the start and end of each group.

    Args:
        L (list): A list of numbers.

    Yields:
        tuple: A tuple representing the start and end of each group of continuous numbers.
    """
    # via https://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
    if len(L) < 1:
        return
    first = last = L[0]
    for n in L[1:]:
        if n - 1 == last:  # Part of the group, bump the end
            last = n
        else:  # Not part of the group, yield current group and start a new
            yield first, last
            first = last = n
    yield first, last  # Yield the last group


def contig_analysis(title, seq):
    """
    Perform analysis on a DNA sequence to identify masked regions and gaps.

    This function examines a DNA sequence to find regions where nucleotides are masked (lowercase) and where gaps (N/n) occur. It uses the `list2groups` function to group consecutive indices of masked regions and gaps.

    Args:
        title (str): The title of the DNA sequence.
        seq (str): The DNA sequence to analyze.

    Returns:
        tuple: A tuple containing:
            - str: The title of the sequence.
            - list: Grouped indices of masked regions.
            - list: Grouped indices of gaps.
    """
    masked = []
    gaps = []
    for i, nuc in enumerate(seq):
        if nuc.islower():
            masked.append(i)
        if nuc in ["N", "n"]:
            gaps.append(i)
    mask_grouped = list(list2groups(masked))
    gaps_grouped = list(list2groups(gaps))
    return (title, mask_grouped, gaps_grouped)


def analyzeAssembly(
    input,
    masked_output,
    gaps_output,
    header_max=16,
    split=False,
    cpus=1,
):
    """
    Analyze the assembly of a genome to identify masked regions and gaps.

    This function processes the genome assembly to detect masked regions and gaps. It utilizes multiprocessing to analyze each contig efficiently and provides statistics on the assembly quality.

    Args:
        input (str): Path to the input genome assembly file.
        masked_output (str): Path to save the softmasked regions output file.
        gaps_output (str): Path to save the assembly gaps output file.
        header_max (int, optional): Maximum length of contig headers. Defaults to 16.
        split (bool, optional): Whether to split contigs into separate files. Defaults to False.
        cpus (int, optional): Number of CPUs to use for multiprocessing. Defaults to 1.

    Returns:
        tuple: A tuple containing:
            - dict: Statistics including number of contigs, assembly size, softmasked percentage, gaps percentage, N50, N90, L50, L90, and average contig length.
            - list: Contig names with headers exceeding the maximum length.
            - list: Errors encountered during analysis.
            - list: Paths to split contig files if split is enabled.
    """
    # local mp functions
    results = []

    def _update(res):
        results.append(res)

    def _handle_error(error):
        print(f"Error: {error}", flush=True)

    # acceptable characters
    IUPAC = {"A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"}

    # load genome into dictionary, build index new just in case
    if os.path.isfile(input + ".fxi"):
        os.remove(input + ".fxi")
    fa = pyfastx.Fasta(input, uppercase=False)

    # get number of each contig and other stats
    bad_names = []
    errors = []
    split_contigs = []
    for base, num in fa.composition.items():
        if base.upper() not in IUPAC:
            errors.append((base, num))

    # now loop through each contig (slow)
    # setup pool
    p = multiprocessing.Pool(cpus)
    for f in fa.keys():
        seq = fa[f]
        if len(seq.name) > header_max:
            bad_names.append(seq.name)
        p.apply_async(
            contig_analysis,
            args=(seq.name, seq.seq),
            callback=_update,
            error_callback=_handle_error,
        )
        if split:
            outname = os.path.join(split, f"{seq.name}.fasta")
            split_contigs.append(outname)
            with open(outname, "w") as outfile:
                outfile.write(f">{seq.name}\n{softwrap(seq.seq)}\n")
    p.close()
    p.join()

    # combine results when finished
    maskLen = 0
    gapLen = 0
    softmasked = {}
    asmgaps = {}
    for r in results:
        title, masked, gaps = r
        maskLen += sum([x[1] - x[0] for x in masked])
        softmasked[title] = masked
        gapLen += sum([y[1] - y[0] for y in gaps])
        asmgaps[title] = gaps

    # calc combined
    percentMask = maskLen / float(fa.size)
    percentGaps = gapLen / float(fa.size)
    if masked_output:
        # write softmasked bed file
        counter = 1
        with open(masked_output, "w") as outfile:
            for k, v in natsorted(list(softmasked.items())):
                for item in v:
                    if len(item) == 2:
                        outfile.write(f"{k}\t{item[0]}\t{item[1]}\tRepeat_{counter}\n")
                        counter += 1
    if gaps_output:
        # write assembly gaps file
        counter2 = 1
        with open(gaps_output, "w") as gapout:
            for k, v in natsorted(list(asmgaps.items())):
                for item in v:
                    if len(item) == 2:
                        gapout.write(
                            f"{k}\t{item[0]}\t{item[1]}\tassembly-gap_{counter2}\n"
                        )
                        counter2 += 1

    # build stats
    n50, l50 = fa.nl(50)
    n90, l90 = fa.nl(90)
    stats = {
        "n_contigs": len(fa),
        "size": fa.size,
        "softmasked": "{:.2%}".format(percentMask),
        "gaps": "{:.2%}".format(percentGaps),
        "n50": n50,
        "n90": n90,
        "l50": l50,
        "l90": l90,
        "avg_length": int(round(fa.mean)),
    }

    return stats, bad_names, errors, split_contigs


def analyzeAssemblySimple(input, header_max=16):
    """
    Analyze the assembly for basic statistics and quality control.

    This function evaluates a genome assembly file to gather basic statistics and identify potential issues. It checks for contig headers exceeding a specified length and identifies any non-standard nucleotide bases.

    Args:
        input (str): Path to the input assembly file.
        header_max (int, optional): Maximum length of contig headers. Defaults to 16.

    Returns:
        tuple: A tuple containing:
            - dict: Statistics including the number of contigs, total size, N50, N90, L50, L90, and average length.
            - list: Contig names exceeding the header length limit.
            - list: Base errors with counts.
    """
    # set character string
    IUPAC = {"A", "C", "G", "T", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"}

    # load genome into dictionary, build index new just in case
    if os.path.isfile(input + ".fxi"):
        os.remove(input + ".fxi")
    fa = pyfastx.Fasta(input, uppercase=False)

    # get number of each contig and other stats
    bad_names = []
    errors = []
    for f in fa.keys():
        seq = fa[f]
        if len(seq.name) > header_max:
            bad_names.append(seq.name)
    for base, num in fa.composition.items():
        if base.upper() not in IUPAC:
            errors.append((base, num))

    # build stats
    n50, l50 = fa.nl(50)
    n90, l90 = fa.nl(90)
    stats = {
        "n_contigs": len(fa),
        "size": fa.size,
        "n50": n50,
        "n90": n90,
        "l50": l50,
        "l90": l90,
        "avg_length": int(round(fa.mean)),
    }

    return stats, bad_names, errors
