import pyhmmer
from .config import env
from .utilities import checkfile
import os
from natsort import natsorted
import subprocess
import json
import json_repair
import uuid
import re


def pyhmmer_version():
    return pyhmmer.__version__


def hmmer_search(hmmfile, sequences, cpus=0, bit_cutoffs=None, evalue=10.0):
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
                            "accession": None
                            if top_hits.query.accession is None
                            else top_hits.query.accession.decode(),
                            "description": None
                            if top_hits.query.description is None
                            else top_hits.query.description.decode(),
                            "seq_length": hit.length,
                            "hmm_length": s_domains[0]["hmm_length"],
                            "hmm_aln_length": s_domains[0]["hmm_aln"],
                            "hmm_coverage": s_domains[0]["hmm_aln"]
                            / s_domains[0]["hmm_length"],
                            "bitscore": hit.score,
                            "evalue": hit.evalue,
                            "domains": s_domains,
                        }
                    )
    return results


def hmmer_scan(hmmfile, sequences, cpus=0, bit_cutoffs=None, evalue=10.0):
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
                            "gene": None
                            if top_hits.query.description is None
                            else top_hits.query.description.decode(),
                            "name": hit.name.decode(),
                            "accession": None
                            if hit.accession is None
                            else hit.accession.decode(),
                            "description": None
                            if hit.description is None
                            else hit.description.decode(),
                            "seq_length": len(top_hits.query.sequence),
                            "hmm_length": s_domains[0]["hmm_length"],
                            "hmm_aln_length": s_domains[0]["hmm_aln"],
                            "hmm_coverage": s_domains[0]["hmm_aln"]
                            / s_domains[0]["hmm_length"],
                            "bitscore": hit.score,
                            "evalue": hit.evalue,
                            "domains": s_domains,
                        }
                    )
    return results


def digitize_sequences(seqs):
    alphabet = pyhmmer.easel.Alphabet.amino()
    sequences = []
    with pyhmmer.easel.SequenceFile(seqs, digital=True, alphabet=alphabet) as seq_file:
        sequences = list(seq_file)
    return sequences


def pfam_search(sequences, cpus=0):
    pfam_hmms = os.path.join(env.get("FUNANNOTATE2_DB"), "Pfam-A.hmm.h3m")
    results = None
    if checkfile(pfam_hmms):
        results = hmmer_search(pfam_hmms, sequences, cpus=cpus, bit_cutoffs="gathering")
    return results


def pfam2tsv(results, output, annots):
    # write the f2 tsv hmm output file as well as the "result"
    with open(output, "w") as outfile:
        json.dump(results, outfile, indent=2)
    with open(annots, "w") as annot:
        for result in natsorted(results, key=lambda x: x["id"]):
            # output the family
            annot.write(f"{result['id']}\tdb_xref\tPFAM:{result['accession']}\n")


def dbcan_search(sequences, cpus=0, evalue=1e-15):
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
            annot.write(f"{result['id']}\tnote\tCAZy:{hit}\n")


def diamond_blast(db, query, cpus=1, evalue=10.0, max_target_seqs=1, tmpdir="/tmp"):
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


def swissprot_blast(
    query, evalue=1e-5, cpus=1, min_pident=60, min_cov=60, max_target_seqs=1
):
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
    with open(output, "w") as outfile:
        json.dump(results, outfile, indent=2)
    with open(annots, "w") as annot:
        for result in results:
            annot.write(
                f"{result['query']}\tdb_xref\tUniProtKB/Swiss-Prot:{result['accession']}\n"
            )


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
