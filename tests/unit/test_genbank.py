"""
Unit tests for the genbank module.
"""
import os
import tempfile
from gfftk.genbank import (
    gff2tbl,
    parse_locus,
    parse_features,
    parse_dbxref,
)


class TestGenbank:
    """Tests for GenBank functions."""

    def test_parse_locus(self):
        """Test parsing GenBank locus information."""
        # Sample locus line
        locus_line = "LOCUS       NC_123456             12345 bp    DNA     linear   CON 01-JAN-2023"

        # Parse the locus line
        result = parse_locus(locus_line)

        # Check the result
        assert result["name"] == "NC_123456"
        assert result["length"] == 12345
        assert result["mol_type"] == "DNA"
        assert result["topology"] == "linear"
        assert result["division"] == "CON"
        assert result["date"] == "01-JAN-2023"

    def test_parse_features_basic(self):
        """Test parsing GenBank features."""
        # Sample feature block
        feature_block = [
            "FEATURES             Location/Qualifiers",
            "     source          1..12345",
            '                     /organism="Test organism"',
            '                     /mol_type="genomic DNA"',
            "     gene            1000..2000",
            '                     /gene="test_gene"',
            '                     /locus_tag="TEST_001"',
            "     CDS             1000..2000",
            '                     /gene="test_gene"',
            '                     /locus_tag="TEST_001"',
            '                     /product="Test protein"',
            '                     /protein_id="XP_123456.1"',
        ]

        # Parse the feature block
        result = parse_features(feature_block)

        # Check the result
        assert len(result) > 0
        assert any(feature["type"] == "gene" for feature in result)
        assert any(feature["type"] == "CDS" for feature in result)

        # Check gene details
        gene = next(feature for feature in result if feature["type"] == "gene")
        assert gene["location"] == "1000..2000"
        assert gene["qualifiers"]["gene"] == ["test_gene"]
        assert gene["qualifiers"]["locus_tag"] == ["TEST_001"]

        # Check CDS details
        cds = next(feature for feature in result if feature["type"] == "CDS")
        assert cds["location"] == "1000..2000"
        assert cds["qualifiers"]["product"] == ["Test protein"]
        assert cds["qualifiers"]["protein_id"] == ["XP_123456.1"]

    def test_parse_dbxref(self):
        """Test parsing database cross-references."""
        # Sample dbxref string
        dbxref_str = "GeneID:12345,NCBI_GP:XP_123456.1"

        # Parse the dbxref string
        result = parse_dbxref(dbxref_str)

        # Check the result
        assert len(result) == 2
        assert result[0] == "GeneID:12345"
        assert result[1] == "NCBI_GP:XP_123456.1"

        # Test with empty string
        assert parse_dbxref("") == []

        # Test with None
        assert parse_dbxref(None) == []

    def test_gff2tbl_basic(self):
        """Test basic GFF to TBL conversion."""
        # Create a temporary GFF file
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".gff3"
        ) as temp_gff:
            temp_gff.write("##gff-version 3\n")
            temp_gff.write(
                "contig1\tprediction\tgene\t1\t1000\t.\t+\t.\tID=gene1;Name=test_gene\n"
            )
            temp_gff.write(
                "contig1\tprediction\tmRNA\t1\t1000\t.\t+\t.\tID=mRNA1;Parent=gene1;Name=test_mrna\n"
            )
            temp_gff.write(
                "contig1\tprediction\texon\t1\t500\t.\t+\t.\tID=exon1;Parent=mRNA1\n"
            )
            temp_gff.write(
                "contig1\tprediction\texon\t600\t1000\t.\t+\t.\tID=exon2;Parent=mRNA1\n"
            )
            temp_gff.write(
                "contig1\tprediction\tCDS\t1\t500\t.\t+\t0\tID=cds1;Parent=mRNA1\n"
            )
            temp_gff.write(
                "contig1\tprediction\tCDS\t600\t1000\t.\t+\t0\tID=cds2;Parent=mRNA1\n"
            )
            temp_gff_name = temp_gff.name

        # Create a temporary FASTA file
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".fasta"
        ) as temp_fasta:
            temp_fasta.write(">contig1\n")
            temp_fasta.write("A" * 1000 + "\n")
            temp_fasta_name = temp_fasta.name

        # Create a temporary output file
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".tbl"
        ) as temp_out:
            temp_out_name = temp_out.name

        try:
            # Convert GFF to TBL
            gff2tbl(temp_gff_name, temp_fasta_name, temp_out_name)

            # Check that the output file exists
            assert os.path.exists(temp_out_name)

            # Check the content of the output file
            with open(temp_out_name, "r") as f:
                content = f.read()

            # Basic checks on the content
            assert ">Feature contig1" in content
            assert "1\t1000\tgene" in content.replace(" ", "")
            assert "gene\ttest_gene" in content
            assert "1\t500" in content
            assert "600\t1000" in content
        finally:
            # Clean up
            for filename in [temp_gff_name, temp_fasta_name, temp_out_name]:
                if os.path.exists(filename):
                    os.unlink(filename)
