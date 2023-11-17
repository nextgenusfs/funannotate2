import sys
import datetime
import random
import warnings
from Bio import BiopythonWarning
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from gfftk.gff import gff2dict
from gfftk.fasta import fasta2dict


def gff3_to_gb(fasta, gff3, gb, offset=600, train_test=False, log=sys.stderr.write):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=BiopythonWarning)

        # write genbank output for augustus training
        # the augustus script is not accurate, results in bad models...
        Seqs = fasta2dict(fasta)
        if isinstance(gff3, str):
            Genes = gff2dict(gff3, fasta, logger=log)
        elif isinstance(gff3, dict):
            Genes = gff3
        # augustus format is non-standard, its each gene as separate record with +/- offset
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
                if v["strand"] == "+":
                    strand = 1
                else:
                    strand = -1
                gene_feat = SeqFeature(
                    FeatureLocation(
                        start=min(v["location"]) - start_pos,
                        end=max(v["location"]) - start_pos + 1,
                        strand=strand,
                    ),
                    type="gene",
                    qualifiers={"locus_tag": k},
                )
                if len(v["mRNA"][0]) > 1:
                    mrna_list = [
                        FeatureLocation(
                            x[0] - start_pos,
                            x[1] - start_pos + 1,
                            strand=strand,
                        )
                        for x in v["mRNA"][0]
                    ]
                    mrna_feat = SeqFeature(
                        CompoundLocation(mrna_list),
                        type="mRNA",
                        qualifiers={"gene": v["ids"][0], "locus_tag": k},
                    )
                else:
                    mrna_feat = SeqFeature(
                        FeatureLocation(
                            start=v["mRNA"][0][0][0] - start_pos,
                            end=v["mRNA"][0][0][1] - start_pos + 1,
                            strand=strand,
                        ),
                        type="mRNA",
                        qualifiers={"gene": v["ids"][0], "locus_tag": k},
                    )
                if len(v["CDS"][0]) > 1:
                    cds_list = [
                        FeatureLocation(
                            x[0] - start_pos, x[1] - start_pos + 1, strand=strand
                        )
                        for x in v["CDS"][0]
                    ]
                    cds_feat = SeqFeature(
                        CompoundLocation(cds_list),
                        type="CDS",
                        qualifiers={
                            "gene": v["ids"][0],
                            "locus_tag": k,
                            "translation": v["protein"][0],
                            "codon_start": v["codon_start"][0],
                        },
                    )
                else:
                    cds_feat = SeqFeature(
                        FeatureLocation(
                            start=v["CDS"][0][0][0] - start_pos,
                            end=v["CDS"][0][0][1] - start_pos + 1,
                            strand=strand,
                        ),
                        type="CDS",
                        qualifiers={
                            "gene": v["ids"][0],
                            "locus_tag": k,
                            "translation": v["protein"][0],
                            "codon_start": v["codon_start"][0],
                        },
                    )
                record = SeqRecord(
                    Seq(sequence),
                    id=f"{v['contig']}_{start_pos}-{end_pos}",
                    annotations={
                        "molecule_type": "DNA",
                        "date": datetime.datetime.now().strftime("%d-%b-%Y"),
                    },
                )
                record.features.append(
                    SeqFeature(FeatureLocation(start=0, end=len(record)), type="source")
                )
                record.features.append(gene_feat)
                record.features.append(mrna_feat)
                record.features.append(cds_feat)
                records.append(record)
        if train_test:
            test_out = f"{gb}.test"
            train_out = f"{gb}.train"
            random.shuffle(records)
            n_recs = int((len(records) + 1) * 0.20)
            if n_recs > 200:
                n_recs = 200
            train_data = records[n_recs:]
            test_data = records[:n_recs]
            with open(test_out, "w") as testout:
                SeqIO.write(test_data, testout, "genbank")
            with open(train_out, "w") as trainout:
                SeqIO.write(train_data, trainout, "genbank")
        with open(gb, "w") as outfile:
            SeqIO.write(records, outfile, "genbank")
        if train_test:
            return len(records), n_recs
        else:
            return len(records), len(records)
