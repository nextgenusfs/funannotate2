import datetime
import random
import sys

import gb_io
from gfftk.fasta import fasta2dict
from gfftk.gff import gff2dict


def gff3_to_gbio(
    fasta, gff3, gb, lowercase=False, offset=600, train_test=False, log=sys.stderr.write
):
    """
    Convert GFF3 annotation to GenBank format.

    This function reads a FASTA file and a GFF3 file (or dictionary), extracts gene features, and converts them into GenBank format. It supports optional lowercase conversion of sequences and splitting data into training and testing sets.

    Args:
        fasta (str): Path to the FASTA file.
        gff3 (str or dict): Path to the GFF3 file or a GFF3 dictionary.
        gb (str): Path to the output GenBank file.
        lowercase (bool, optional): Whether to convert the sequence to lowercase. Defaults to False.
        offset (int, optional): Offset value for sequence positioning. Defaults to 600.
        train_test (bool, optional): Whether to split the data into training and testing sets. Defaults to False.
        log (function, optional): Logger function. Defaults to sys.stderr.write.

    Returns:
        tuple: Number of total records and number of records in the test set if train_test is True, otherwise returns the number of total records twice.
    """

    def _fetch_coords(v, i=0, start_pos=0, feature="gene"):
        if feature == "gene":
            coords = gb_io.Range(
                min(v["location"]) - start_pos,
                max(v["location"]) - start_pos + 1,
            )
            if v["strand"] == "-":
                return gb_io.Complement(coords)
            else:
                return coords
        else:
            if feature == "CDS":
                if len(v["CDS"][i]) > 1:
                    coords = gb_io.Join(
                        [
                            gb_io.Range(x[0] - start_pos, x[1] - start_pos + 1)
                            for x in sorted(v["CDS"][i])
                        ]
                    )
                else:
                    coords = gb_io.Range(
                        v["CDS"][i][0][0] - start_pos,
                        v["CDS"][i][0][1] - start_pos + 1,
                    )
                if v["strand"] == "-":
                    return gb_io.Complement(coords)
                else:
                    return coords
            else:
                if len(v["mRNA"][i]) > 1:
                    coords = gb_io.Join(
                        [
                            gb_io.Range(x[0] - start_pos, x[1] - start_pos + 1)
                            for x in sorted(v["mRNA"][i])
                        ]
                    )
                else:
                    coords = gb_io.Range(
                        v["mRNA"][i][0][0] - start_pos,
                        v["mRNA"][i][0][1] - start_pos + 1,
                    )
                if v["strand"] == "-":
                    return gb_io.Complement(coords)
                else:
                    return coords

    Seqs = fasta2dict(fasta, full_header=True)
    if isinstance(gff3, str):
        Genes = gff2dict(gff3, fasta, logger=log)
    elif isinstance(gff3, dict):
        Genes = gff3
    records = []
    for k, v in Genes.items():
        loff = offset
        roff = offset
        if (
            "mRNA" in v["type"]
            and len(v["CDS"][0]) > 0
            and v["protein"][0].startswith("M")
            and v["protein"][0].endswith("*")
        ):
            start_pos = min(v["location"]) - loff
            if start_pos <= 0:
                start_pos = 1
                start_idx = 0
                loff = min(v["location"])
            else:
                start_idx = start_pos - 1
            end_pos = max(v["location"]) + roff
            if end_pos > len(Seqs[v["contig"]]):
                end_pos = len(Seqs[v["contig"]])
            sequence = Seqs[v["contig"]][start_idx:end_pos]
            # now collect features
            rec_features = [
                gb_io.Feature(
                    "source",
                    gb_io.Range(0, len(sequence)),
                    qualifiers=[
                        gb_io.Qualifier("mol_type", value="genomic DNA"),
                    ],
                ),
            ]
            # gene feature
            rec_features.append(
                gb_io.Feature(
                    "gene",
                    _fetch_coords(v, start_pos=start_pos, feature="gene"),
                    qualifiers=[gb_io.Qualifier("locus_tag", value=k)],
                )
            )
            # transcript feature
            rec_features.append(
                gb_io.Feature(
                    "mRNA",
                    _fetch_coords(v, start_pos=start_pos, feature="mRNA"),
                    qualifiers=[
                        gb_io.Qualifier("gene", value=v["ids"][0]),
                        gb_io.Qualifier("locus_tag", value=k),
                    ],
                )
            )
            # coding sequence feature
            rec_features.append(
                gb_io.Feature(
                    "CDS",
                    _fetch_coords(v, start_pos=start_pos, feature="CDS"),
                    qualifiers=[
                        gb_io.Qualifier("gene", value=v["ids"][0]),
                        gb_io.Qualifier("locus_tag", value=k),
                        gb_io.Qualifier("translation", value=v["protein"][0]),
                        gb_io.Qualifier("codon_start", value=str(v["codon_start"][0])),
                    ],
                )
            )
            if lowercase:
                sequence = sequence.lower()
            records.append(
                gb_io.Record(
                    sequence=str.encode(sequence),
                    circular=False,
                    molecule_type="DNA",
                    keywords=".",
                    definition=".",
                    accession=f"{v['contig']}_{start_pos}-{end_pos}",
                    name=f"{v['contig']}_{start_pos}-{end_pos}",
                    version=f"{v['contig']}_{start_pos}-{end_pos}",
                    length=len(sequence),
                    date=datetime.date.today(),
                    features=rec_features,
                )
            )
    if train_test:
        test_out = f"{gb}.test"
        train_out = f"{gb}.train"
        random.shuffle(records)
        n_recs = int((len(records) + 1) * 0.20)
        if n_recs > 200:
            n_recs = 200
        train_data = records[n_recs:]
        test_data = records[:n_recs]
        with open(test_out, "wb") as testout:
            gb_io.dump(test_data, testout, escape_locus=False, truncate_locus=False)
        with open(train_out, "wb") as trainout:
            gb_io.dump(train_data, trainout, escape_locus=False, truncate_locus=False)

    with open(gb, "wb") as outfile:
        gb_io.dump(records, outfile, escape_locus=False, truncate_locus=False)
    if train_test:
        return len(records), n_recs
    else:
        return len(records), len(records)
