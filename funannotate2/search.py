import json
import os
import re
import subprocess
import sys
import uuid

import json_repair
import pyhmmer
from buscolite.busco import runbusco
from natsort import natsorted

from .config import env
from .utilities import checkfile


def pyhmmer_version():
    return pyhmmer.__version__


def hmmer_search(hmmfile, sequences, cpus=0, bit_cutoffs=None, evalue=10.0):
    """
    Perform an HMMER search using the provided HMM file and sequences.

    This function utilizes the `pyhmmer` library to search for sequence alignments against a given HMM file. It processes the results to extract relevant information about each hit, including domain details and scoring metrics.

    Args:
        hmmfile (str): Path to the HMM file.
        sequences (list): List of sequences to search against.
        cpus (int, optional): Number of CPUs to use for the search. Defaults to 0.
        bit_cutoffs (list, optional): List of bit score cutoffs for reporting hits. Defaults to None.
        evalue (float, optional): E-value threshold for reporting hits. Defaults to 10.0.

    Returns:
        list: A list of dictionaries containing information about the hits found, including:
            - id (str): Identifier of the hit.
            - gene (str): Gene description of the hit.
            - name (str): Name of the query sequence.
            - accession (str or None): Accession of the query sequence.
            - description (str or None): Description of the query sequence.
            - seq_length (int): Length of the sequence.
            - hmm_length (int): Length of the HMM.
            - hmm_aln_length (int): Alignment length of the HMM.
            - hmm_coverage (float): Coverage of the HMM alignment.
            - bitscore (float): Bit score of the hit.
            - evalue (float): E-value of the hit.
            - domains (list): List of domain information dictionaries.
    """
    results = []
    with pyhmmer.plan7.HMMFile(hmmfile) as hfile:
        for top_hits in pyhmmer.hmmsearch(
            hfile, sequences, cpus=cpus, bit_cutoffs=bit_cutoffs, E=evalue
        ):
            for hit in top_hits:
                if hit.included:
                    domains = []
                    for h in hit.domains:
                        domains.append(
                            {
                                "hmm_from": h.alignment.hmm_from,
                                "hmm_to": h.alignment.hmm_to,
                                "hmm_length": h.alignment.hmm_length,
                                "hmm_aln": h.alignment.hmm_to - h.alignment.hmm_from,
                                "env_from": h.env_from,
                                "env_to": h.env_to,
                                "env_length": h.env_to - h.env_from,
                                "score": h.score,
                            }
                        )
                    # sort by highest scoring domain hit
                    s_domains = sorted(domains, key=lambda x: x["score"], reverse=True)
                    results.append(
                        {
                            "id": hit.name.decode(),
                            "gene": hit.description.decode(),
                            "name": top_hits.query.name.decode(),
                            "accession": (
                                None
                                if top_hits.query.accession is None
                                else top_hits.query.accession.decode()
                            ),
                            "description": (
                                None
                                if top_hits.query.description is None
                                else top_hits.query.description.decode()
                            ),
                            "seq_length": hit.length,
                            "hmm_length": s_domains[0]["hmm_length"],
                            "hmm_aln_length": s_domains[0]["hmm_aln"],
                            "hmm_coverage": s_domains[0]["hmm_aln"] / s_domains[0]["hmm_length"],
                            "bitscore": hit.score,
                            "evalue": hit.evalue,
                            "domains": s_domains,
                        }
                    )
    return results


def hmmer_scan(hmmfile, sequences, cpus=0, bit_cutoffs=None, evalue=10.0):
    """
    Perform HMMER scan on sequences using the provided HMM file.

    This function utilizes the `pyhmmer` library to scan sequences against a specified HMM file. It processes the results to extract relevant information about each hit, including domain details and scoring metrics.

    Args:
        hmmfile (str): Path to the HMM file.
        sequences (list): List of sequences to scan.
        cpus (int, optional): Number of CPUs to use. Defaults to 0.
        bit_cutoffs (float, optional): Bit score cutoff for reporting hits. Defaults to None.
        evalue (float, optional): E-value threshold for reporting hits. Defaults to 10.0.

    Returns:
        list: A list of dictionaries containing information about the hits, including:
            - id (str): Identifier of the query sequence.
            - gene (str or None): Gene description of the query sequence.
            - name (str): Name of the hit.
            - accession (str or None): Accession of the hit.
            - description (str or None): Description of the hit.
            - seq_length (int): Length of the query sequence.
            - hmm_length (int): Length of the HMM.
            - hmm_aln_length (int): Alignment length of the HMM.
            - hmm_coverage (float): Coverage of the HMM alignment.
            - bitscore (float): Bit score of the hit.
            - evalue (float): E-value of the hit.
            - domains (list): List of domain information dictionaries.
    """
    results = []
    with pyhmmer.plan7.HMMFile(hmmfile) as hfile:
        for top_hits in pyhmmer.hmmscan(
            sequences, hfile, cpus=cpus, bit_cutoffs=bit_cutoffs, E=evalue
        ):
            for hit in top_hits:
                if hit.included:
                    domains = []
                    for h in hit.domains:
                        domains.append(
                            {
                                "hmm_from": h.alignment.hmm_from,
                                "hmm_to": h.alignment.hmm_to,
                                "hmm_length": h.alignment.hmm_length,
                                "hmm_aln": h.alignment.hmm_to - h.alignment.hmm_from,
                                "env_from": h.env_from,
                                "env_to": h.env_to,
                                "env_length": h.env_to - h.env_from,
                                "score": h.score,
                            }
                        )
                    # sort by highest scoring domain hit
                    s_domains = sorted(domains, key=lambda x: x["score"], reverse=True)
                    results.append(
                        {
                            "id": top_hits.query.name.decode(),
                            "gene": (
                                None
                                if top_hits.query.description is None
                                else top_hits.query.description.decode()
                            ),
                            "name": hit.name.decode(),
                            "accession": (
                                None if hit.accession is None else hit.accession.decode()
                            ),
                            "description": (
                                None if hit.description is None else hit.description.decode()
                            ),
                            "seq_length": len(top_hits.query.sequence),
                            "hmm_length": s_domains[0]["hmm_length"],
                            "hmm_aln_length": s_domains[0]["hmm_aln"],
                            "hmm_coverage": s_domains[0]["hmm_aln"] / s_domains[0]["hmm_length"],
                            "bitscore": hit.score,
                            "evalue": hit.evalue,
                            "domains": s_domains,
                        }
                    )
    return results


def digitize_sequences(seqs):
    """
    Digitize a list of sequences using the pyhmmer library.

    This function reads a sequence file and converts the sequences into a digital format using the `pyhmmer.easel` module. It utilizes an amino acid alphabet for the conversion.

    Args:
        seqs (str): Path to the sequence file to be digitized.

    Returns:
        list: A list of digitized sequences.
    """
    alphabet = pyhmmer.easel.Alphabet.amino()
    sequences = []
    with pyhmmer.easel.SequenceFile(seqs, digital=True, alphabet=alphabet) as seq_file:
        sequences = list(seq_file)
    return sequences


def pfam_search(sequences, cpus=0):
    """
    Perform a Pfam search using the Pfam-A HMM database and the provided sequences.

    This function checks for the existence of the Pfam-A HMM database file and uses the `hmmer_search` function to search for sequence alignments. It utilizes the `checkfile` utility to verify the file's presence and size.

    Args:
        sequences (list): List of sequences to search against.
        cpus (int, optional): Number of CPUs to use for the search. Defaults to 0.

    Returns:
        list: A list of dictionaries containing information about the hits found using the Pfam search, or None if the database file is not found.
    """
    pfam_hmms = os.path.join(env.get("FUNANNOTATE2_DB"), "Pfam-A.hmm.h3m")
    results = None
    if checkfile(pfam_hmms):
        results = hmmer_search(pfam_hmms, sequences, cpus=cpus, bit_cutoffs="gathering")
    return results


def pfam2tsv(results, output, annots):
    """
    Write PFAM results to a TSV file and generate annotations.

    This function processes PFAM search results, writes them to a TSV file, and generates annotations in a separate file. It uses the `add2dict` utility to organize annotations in a dictionary format.

    Args:
        results (list): List of PFAM results.
        output (str): Path to the output TSV file.
        annots (str): Path to the annotations file.

    Returns:
        dict: A dictionary containing annotations, where keys are gene IDs and values are lists of annotation details.
    """
    a = {}
    # write the f2 tsv hmm output file as well as the "result"
    with open(output, "w") as outfile:
        json.dump(results, outfile, indent=2)
    with open(annots, "w") as annot:
        for result in natsorted(results, key=lambda x: x["id"]):
            # output the family
            annot.write(f"{result['id']}\tdb_xref\tPFAM:{result['accession']}\n")
            a = add2dict(a, result["id"], "db_xref", f"PFAM:{result['accession']}")
    return a


def dbcan_search(sequences, cpus=0, evalue=1e-15):
    """
    Perform a search using dbCAN HMM profiles on the provided sequences.

    This function utilizes the `hmmer_scan` function to search sequences against the dbCAN HMM database. It filters the results based on HMM coverage and organizes them by subfamily, selecting the best hit for each subfamily.

    Args:
        sequences (list): List of sequences to search.
        cpus (int, optional): Number of CPUs to use. Defaults to 0.
        evalue (float, optional): E-value threshold for reporting hits. Defaults to 1e-15.

    Returns:
        list: A list of dictionaries containing information about the best hits for each sequence, filtered by HMM coverage.
    """
    dbcan_hmms = os.path.join(env.get("FUNANNOTATE2_DB"), "dbCAN.hmm.h3m")
    results = None
    if checkfile(dbcan_hmms):
        raw_results = hmmer_scan(dbcan_hmms, sequences, cpus=cpus, evalue=evalue)
        results = []
        filt = {}
        # dbcan should be filtered by hmm length, greater than 0.35
        for res in raw_results:
            if res["hmm_coverage"] < 0.35:
                continue
            family = re.match(r"^[A-Z]+", res["name"]).group(0)
            if res["id"] not in filt:
                filt[res["id"]] = {family: [res]}
            else:
                if family not in filt[res["id"]]:
                    filt[res["id"]][family] = [res]
                else:
                    filt[res["id"]][family].append(res)
        # now get best subfamily hit for each
        for k, v in filt.items():
            for f, resL in v.items():
                sl = sorted(resL, key=lambda x: x["evalue"])
                results.append(sl[0])
    return results


def dbcan2tsv(results, output, annots):
    """
    Write dbCAN results to a TSV file and generate annotations.

    This function processes dbCAN search results, writes them to a TSV file, and generates annotations in a separate file. It uses the `add2dict` utility to organize annotations in a dictionary format.

    Args:
        results (list): List of dbCAN results.
        output (str): Path to the output TSV file.
        annots (str): Path to the annotations file.

    Returns:
        dict: A dictionary containing annotations, where keys are gene IDs and values are lists of annotation details.
    """
    a = {}
    # write the f2 tsv hmm output file as well as the "result"
    with open(output, "w") as outfile:
        json.dump(results, outfile, indent=2)
    with open(annots, "w") as annot:
        for result in natsorted(results, key=lambda x: x["id"]):
            # output the family
            if "_" in result["name"]:
                hit = result["name"].split("_")[0]
            else:
                hit = result["name"]
            hit = hit.replace(".hmm", "")
            annot.write(f"{result['id']}\tnote\tCAZy:{hit}\n")
            a = add2dict(a, result["id"], "note", f"CAZy:{hit}")
    return a


def diamond_blast(db, query, cpus=1, evalue=10.0, max_target_seqs=1, tmpdir="/tmp"):
    """
    Run a protein BLAST using Diamond.

    This function executes a Diamond BLAST search with the specified parameters, processes the output in JSON format, and returns the results. It uses a temporary file to store intermediate results and ensures cleanup after processing.

    Args:
        db (str): Path to the Diamond database.
        query (str): Path to the query file.
        cpus (int, optional): Number of CPUs to use. Defaults to 1.
        evalue (float, optional): E-value threshold. Defaults to 10.0.
        max_target_seqs (int, optional): Maximum number of target sequences to report. Defaults to 1.
        tmpdir (str, optional): Temporary directory to store intermediate files. Defaults to "/tmp".

    Returns:
        list: JSON formatted results of the BLAST search.
    """
    # use diamond as protein blast
    tmpfile = os.path.join(tmpdir, f"diamond_{uuid.uuid4()}")
    cmd = [
        "diamond",
        "blastp",
        "--sensitive",
        "--query",
        query,
        "--threads",
        str(cpus),
        "--db",
        db,
        "--evalue",
        str(evalue),
        "--max-target-seqs",
        str(max_target_seqs),
        "--outfmt",
        "104",
        "qseqid",
        "qtitle",
        "qlen",
        "sseqid",
        "stitle",
        "slen",
        "pident",
        "length",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "--out",
        tmpfile,
    ]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.communicate()
    with open(tmpfile) as infile:
        results = json_repair.load(infile)
    if os.path.isfile(tmpfile):
        os.remove(tmpfile)
    return results


def merops_blast(query, evalue=1e-5, cpus=1, max_target_seqs=1):
    """
    Run a BLAST search against the MEROPS database using Diamond.

    This function performs a BLAST search against the MEROPS database with the specified parameters. It calculates coverage for each result and extracts family information from the search output.

    Args:
        query (str): Path to the query file.
        evalue (float, optional): E-value threshold for reporting hits. Defaults to 1e-5.
        cpus (int, optional): Number of CPUs to use. Defaults to 1.
        max_target_seqs (int, optional): Maximum number of target sequences to report. Defaults to 1.

    Returns:
        list: A list of dictionaries containing information about the hits found, including coverage and family details.
    """
    merops_db = os.path.join(env.get("FUNANNOTATE2_DB"), "merops.dmnd")
    results = None
    if checkfile(merops_db):
        results = []
        raw_results = diamond_blast(
            merops_db, query, cpus=cpus, evalue=evalue, max_target_seqs=max_target_seqs
        )
        for res in raw_results:
            coverage = (res["length"] / res["qlen"]) * 100
            res["cov"] = coverage
            mer_id, mer_family = res["stitle"].split()
            res["family"] = mer_family
            results.append(res)
    return results


def merops2tsv(results, output, annots):
    """
    Write MEROPS results to TSV files and generate annotations.

    This function processes MEROPS search results, writes them to a TSV file, and generates annotations in a separate file. It uses the `add2dict` utility to organize annotations in a dictionary format.

    Args:
        results (list): List of MEROPS search results.
        output (str): Path to the output file for JSON results.
        annots (str): Path to the output file for annotations.

    Returns:
        dict: A dictionary containing annotations, where keys are query sequence IDs and values are lists of annotation details.
    """
    a = {}
    with open(output, "w") as outfile:
        json.dump(results, outfile, indent=2)
    with open(annots, "w") as annot:
        for result in results:
            annot.write(f"{result['qseqid']}\tnote\tMEROPS:{result['sseqid']} {result['family']}\n")
            a = add2dict(
                a,
                result["qseqid"],
                "note",
                f"MEROPS:{result['sseqid']} {result['family']}",
            )
    return a


def swissprot_blast(query, evalue=1e-5, cpus=1, min_pident=60, min_cov=60, max_target_seqs=1):
    """
    Perform a BLAST search against the SwissProt database using Diamond.

    This function uses the Diamond tool to search a query file against the SwissProt database. It filters results based on minimum percentage identity and coverage, and parses the SwissProt headers for additional information.

    Args:
        query (str): Path to the query file.
        evalue (float, optional): E-value threshold for reporting hits. Defaults to 1e-5.
        cpus (int, optional): Number of CPUs to use. Defaults to 1.
        min_pident (int, optional): Minimum percentage identity required for a hit. Defaults to 60.
        min_cov (int, optional): Minimum coverage percentage required for a hit. Defaults to 60.
        max_target_seqs (int, optional): Maximum number of target sequences to report. Defaults to 1.

    Returns:
        list: A list of dictionaries containing information about the filtered hits, including coverage, percentage identity, and alignment details.
    """
    uniprot_db = os.path.join(env.get("FUNANNOTATE2_DB"), "uniprot.dmnd")
    results = None
    if checkfile(uniprot_db):
        results = []
        raw_results = diamond_blast(
            uniprot_db, query, cpus=cpus, evalue=evalue, max_target_seqs=max_target_seqs
        )
        for res in raw_results:
            if res["pident"] >= min_pident:
                coverage = (res["length"] / res["qlen"]) * 100
                if coverage >= min_cov:
                    r = parse_swissprot_headers(res["stitle"])
                    r["cov"] = round(coverage, 2)
                    r["pident"] = res["pident"]
                    r["evalue"] = res["evalue"]
                    r["aln_length"] = res["length"]
                    r["query"] = res["qseqid"]
                    r["gene"] = res["qtitle"].split()[-1]
                    r["qstart"] = res["qstart"]
                    r["qend"] = res["qend"]
                    r["sstart"] = res["sstart"]
                    r["send"] = res["send"]
                    r["slen"] = res["slen"]
                    r["qlen"] = res["qlen"]
                    results.append(r)
    return results


def parse_swissprot_headers(header):
    """
    Parse the header of a SwissProt entry to extract relevant information.

    This function processes a SwissProt header string to extract and organize details such as the description, accession number, species name, and additional data fields. It handles specific formatting and filtering of unwanted words in the description.

    Args:
        header (str): The header string of the SwissProt entry.

    Returns:
        dict: A dictionary containing the parsed information, including:
            - description (str): The cleaned description of the entry.
            - accession (str): The accession number of the entry.
            - sp_name (str): The species name from the entry.
            - Additional fields extracted from the header, such as gene name (GN) and others.
    """
    # first space is the id
    sseqid, rest = header.split(" ", 1)
    _, accession, sp_name = sseqid.split("|")
    desc = re.split(r"\ [A-Z]{2}\=", rest)[0]
    # fix description
    bad_desc_words = ["(Fragment)", "homolog,", "AltName:"]
    tmp_desc = [x for x in desc.split() if x not in bad_desc_words]
    description = " ".join(tmp_desc)
    keyed_data = rest.split(desc)[-1]
    split_data = [x.strip() for x in re.split(r"[A-Z]{2}[=]{1}", keyed_data)[1:]]
    d = {"description": description, "accession": accession, "sp_name": sp_name}
    for i, n in enumerate(re.findall(r"[A-Z]{2}[=]{1}", keyed_data)):
        if n == "GN=":
            raw = split_data[i]
            raw = raw.split(" ")[0].upper().replace("-", "")
            d[n.strip("=")] = raw
        else:
            d[n.rstrip("=")] = split_data[i]
    return d


def swissprot2tsv(results, output, annots):
    """
    Converts Swiss-Prot results to TSV format and generates annotations.

    This function processes a list of Swiss-Prot results, writes the results to a JSON file,
    and generates annotations in a TSV format. It returns a dictionary containing the processed
    results with additional annotations.

    Parameters:
    - results (list): A list of dictionaries containing Swiss-Prot result data.
    - output (str): The file path where the JSON output will be written.
    - annots (str): The file path where the TSV annotations will be written.

    Returns:
    - dict: A dictionary containing the processed results with annotations.
    """
    # return a dictionary of results
    a = {}
    with open(output, "w") as outfile:
        json.dump(results, outfile, indent=2)
    with open(annots, "w") as annot:
        for result in results:
            annot.write(f"{result['query']}\tdb_xref\tUniProtKB/Swiss-Prot:{result['accession']}\n")
            # add db_xref
            a = add2dict(
                a,
                result["query"],
                "db_xref",
                f"UniProtKB/Swiss-Prot:{result['accession']}",
            )
            # try to fix description
            product = swissprot_clean_product(result["description"])
            # see if gene name is valid
            if "GN" in result:
                if swissprot_valid_gene(result["GN"]):
                    a = add2dict(a, result["query"], "name", result["GN"])
                    annot.write(f"{result['query']}\tname\t{result['GN']}\n")
                    a = add2dict(a, result["query"], "product", product)
                    annot.write(f"{result['query']}\tproduct\t{product}\n")
                else:
                    # add as a note
                    a = add2dict(
                        a,
                        result["query"],
                        "note",
                        f"{result['pident']}% identical to {result['sp_name']} {product}",
                    )
                    annot.write(
                        f"{result['query']}\tnote\t{result['pident']}% identical to {result['sp_name']} {product}\n"
                    )
            else:
                # add as a note
                a = add2dict(
                    a,
                    result["query"],
                    "note",
                    f"{result['pident']}% identical to {result['sp_name']} {product}",
                )
                annot.write(
                    f"{result['query']}\tnote\t{result['pident']}% identical to {result['sp_name']} {product}\n"
                )
    return a


def parse_annotations(tsv):
    """
    Parse a three-column annotation file into a dictionary.

    This function reads a TSV file containing annotations, where each line consists of a gene,
    a database identifier, and a value. It processes the file and returns a dictionary where
    each gene is a key, and its associated database and value are stored as entries.

    Parameters:
    - tsv (str): The file path to the annotation file in TSV format.

    Returns:
    - dict: A dictionary with genes as keys and their corresponding database and value entries.
    """
    # parse a three column annotation file into a dictionary
    a = {}
    with open(tsv, "r") as infile:
        for line in infile:
            line = line.rstrip()
            gene, db, value = line.split("\t")
            a = add2dict(a, gene, db, value)
    return a


def add2dict(adict, gene, key, value):
    """
    Add a key-value pair to a nested dictionary structure.

    This function updates a dictionary by adding a value to a list associated with a specific key
    within a nested dictionary. If the gene is not present in the dictionary, a new entry is created.
    If the key already exists for the gene, the value is appended to the list; otherwise, a new list
    is created for the key.

    Parameters:
    - adict (dict): The dictionary to be updated.
    - gene (str): The key for the outer dictionary.
    - key (str): The key for the inner dictionary.
    - value (any): The value to be added to the list associated with the inner key.

    Returns:
    - dict: The updated dictionary with the new key-value pair.
    """
    if gene not in adict:
        adict[gene] = {key: [value]}
    else:
        if key in adict[gene]:
            adict[gene][key].append(value)
        else:
            adict[gene][key] = [value]
    return adict


def swissprot_valid_gene(name):
    if number_present(name) and len(name) > 2 and not morethanXnumbers(name, 3) and "." not in name:
        return True
    else:
        return False


def swissprot_clean_product(description):
    """
    Clean a product description by removing specified unwanted words.

    This function processes a product description string by removing certain predefined words
    such as "(Fragment)", "homolog", "homolog,", and "AltName:". The cleaned description is
    returned as a string without these unwanted words.

    Parameters:
    - description (str): The input product description to be cleaned.

    Returns:
    - str: The cleaned product description with specified words removed.
    """
    # need to do some filtering here of certain words
    bad_words = ["(Fragment)", "homolog", "homolog,", "AltName:"]
    # turn string into array, splitting on spaces
    descript = description.split(" ")
    final_desc = [x for x in descript if x not in bad_words]
    final_desc = " ".join(final_desc)
    return final_desc


def number_present(s):
    return any(i.isdigit() for i in s)


def morethanXnumbers(s, num):
    count = 0
    for i in s:
        if number_present(i):
            count += 1
    if count >= num:
        return True
    else:
        return False


def capfirst(x):
    return x[0].upper() + x[1:]


def busco_search(query, buscodb, cpus=1, logger=sys.stderr):
    """
    Run BUSCO search on a query using a specified BUSCO database.

    This function executes a BUSCO search on the provided query file using the specified BUSCO database.
    It utilizes the `runbusco` function to perform the search in protein mode and returns the results.

    Args:
        query (str): Path to the query file to be analyzed.
        buscodb (str): Path to the BUSCO database to be used for the search.
        cpus (int, optional): Number of CPUs to allocate for the search process. Defaults to 1.
        logger (file object, optional): Logger object for logging messages. Defaults to sys.stderr.

    Returns:
        dict: The results of the BUSCO search, containing information about the identified BUSCOs.
    """
    busco_results, busco_missing, stats, cfg = runbusco(
        query,
        buscodb,
        mode="proteins",
        cpus=cpus,
        logger=logger,
        verbosity=0,
    )
    return busco_results


def busco2tsv(results, buscodb, busco_results, annots):
    """
    Parse BUSCO proteome analysis results and construct functional annotations.

    This function processes the results of a BUSCO analysis by writing the raw results to a file,
    extracting database information, and constructing functional annotations. The annotations are
    written to a specified file in a tab-separated format.

    Parameters:
    - results (dict): The BUSCO analysis results, where each key is a BUSCO ID and the value is a dictionary
      containing details such as 'name', 'hit', 'bitscore', 'evalue', 'domains', 'length', and 'status'.
    - buscodb (str): The path to the BUSCO database directory, which contains links to additional information.
    - busco_results (str): The file path where the raw BUSCO results will be written in JSON format.
    - annots (str): The file path where the constructed functional annotations will be written.

    Returns:
    - dict: A dictionary containing the constructed functional annotations, where each key is a gene hit
      and the value is a dictionary with annotation details.
    """
    # function to parse busco proteome analysis
    # first write the raw data
    with open(busco_results, "w") as outfile:
        json.dump(results, outfile, indent=2)
    # get the DB info and load for constructing the functional annotation
    odb_version = os.path.basename(buscodb)
    links_file = None
    for f in os.listdir(buscodb):
        if f.startswith("links_to_ODB"):
            links_file = os.path.join(buscodb, f)
    busco_data = {}
    if links_file:
        with open(links_file, "r") as infile:
            for line in infile:
                line = line.strip()
                busco_id, description, url = line.split("\t")
                busco_data[busco_id] = {"description": description, "url": url}
    # the results object is a dictionary like this:
    # 16968at33183 {'name': '16968at33183', 'hit': 'FUN2_002340-T1', 'bitscore': 148.97767639160156, 'evalue': 1.1007091320326609e-44, 'domains': [{'hmm_from': 2, 'hmm_to': 51, 'env_from': 2, 'env_to': 60, 'score': 70.0794906616211}, {'hmm_from': 46, 'hmm_to': 97, 'env_from': 54, 'env_to': 132, 'score': 82.8307876586914}], 'length': 133, 'status': 'complete'}
    # so now we want to construct the 3 column
    a = {}
    with open(annots, "w") as annot:
        for k, v in natsorted(results.items(), key=lambda x: (x[1]["hit"], -x[1]["bitscore"])):
            if v["name"] in busco_data:
                defline = f'BUSCO:{k} [{odb_version}] {busco_data.get(v["name"])["description"]}'
            else:
                defline = f"BUSCO:{k} [{odb_version}]"
            if v["hit"] not in a:
                a = add2dict(a, v["hit"], "note", defline)
                annot.write(f"{v['hit']}\tnote\t{defline}\n")
    return a
