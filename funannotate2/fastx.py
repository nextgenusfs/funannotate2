import pyfastx
from natsort import natsorted
import os
import concurrent.futures
from .utilities import runProcessJob
import multiprocessing


def countfasta(fasta_file):
    count = 0
    for seq in pyfastx.Fasta(fasta_file, build_index=False):
        count += 1
    return count


def softwrap(string, every=80):
    lines = []
    for i in range(0, len(string), every):
        lines.append(string[i : i + every])
    return "\n".join(lines)


def fasta2dict(fasta_file):
    fa = {}
    for title, seq in pyfastx.Fasta(fasta_file, build_index=False):
        fa[title] = seq
    return fa


def list2groups(L):
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
    masked = []
    gaps = []
    for i, nuc in enumerate(seq):
        if nuc.islower():
            masked.append(i)
        nuc = nuc.upper()
        if nuc == "N":
            masked.append(i)
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
        maskLen += len(masked)
        softmasked[title] = masked
        gapLen += len(gaps)
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
        if base not in IUPAC:
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
