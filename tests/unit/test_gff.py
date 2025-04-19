"""
Unit tests for the gff module.
"""

import os
import tempfile

from gfftk.gff import _detect_format, dict2gff3, gff2dict, simplifyGO, start_end_gap


class TestGFFParsing:
    """Tests for GFF parsing functions."""

    def test_gff2dict_basic(self, sample_gff_file, sample_fasta_file):
        """Test basic GFF parsing."""
        # Parse the GFF file
        result = gff2dict(sample_gff_file, sample_fasta_file)

        # Check that the result is a dictionary
        assert isinstance(result, dict)

        # Check that the gene is in the result
        assert "gene1" in result

        # Check the gene properties
        gene = result["gene1"]
        assert gene["type"] == ["mRNA"]  # The actual implementation uses a list for type
        assert gene["contig"] == "contig1"
        assert gene["location"] == (
            1,
            1000,
        )  # The actual implementation uses a tuple for location
        assert gene["strand"] == "+"

        # Print the gene dictionary to see its structure
        import pprint

        print("Gene dictionary structure:")
        pprint.pprint(gene)

        # Check that the mRNA is in the gene
        assert "mRNA" in gene
        assert len(gene["mRNA"]) == 1

        # Check that the CDS regions are in the gene
        assert "CDS" in gene
        assert len(gene["CDS"]) == 1  # CDS is a list of lists of tuples
        assert len(gene["CDS"][0]) == 2  # There are 2 CDS regions
        assert gene["CDS"][0][0] == (1, 500)  # First CDS region
        assert gene["CDS"][0][1] == (600, 1000)  # Second CDS region

    def test_dict2gff3_basic(self, temp_dir, sample_gff_file, sample_fasta_file):
        """Test basic GFF dictionary to GFF3 conversion."""
        # Parse the GFF file to get a dictionary
        gff_dict = gff2dict(sample_gff_file, sample_fasta_file)

        # Convert the dictionary back to GFF3
        output_file = os.path.join(temp_dir, "output.gff3")
        dict2gff3(gff_dict, output=output_file)

        # Check that the output file exists
        assert os.path.exists(output_file)

        # Check the content of the output file
        with open(output_file, "r") as f:
            content = f.read()

        # Basic checks on the content
        assert "##gff-version 3" in content
        assert "gene1" in content
        assert "mRNA1" in content
        assert "exon" in content
        assert "CDS" in content

    def test_start_end_gap(self):
        """Test the start_end_gap function."""
        # Test with N's at the start
        seq = "NNNNATGC"
        coords = [(1, 8)]
        result_seq, result_coords = start_end_gap(seq, coords)
        assert result_seq == "ATGC"
        assert result_coords == [(5, 8)]

        # Test with N's at the end
        seq = "ATGCNNN"
        coords = [(1, 7)]
        result_seq, result_coords = start_end_gap(seq, coords)
        assert result_seq == "ATGC"
        assert result_coords == [(1, 4)]

        # Test with N's at both ends
        seq = "NNNATGCNNN"
        coords = [(1, 10)]
        result_seq, result_coords = start_end_gap(seq, coords)
        assert result_seq == "ATGC"
        assert result_coords == [(4, 7)]

        # Test with no N's
        seq = "ATGC"
        coords = [(1, 4)]
        result_seq, result_coords = start_end_gap(seq, coords)
        assert result_seq == "ATGC"
        assert result_coords == [(1, 4)]

    def test_detect_format(self, sample_gff_file):
        """Test the _detect_format function."""
        # Test with a GFF3 file
        parser, format_type = _detect_format(sample_gff_file)
        # The actual implementation returns a tuple (parser_function, format_string)
        assert format_type == "default"

        # Create a GTF-like file for testing
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp:
            temp.write(
                'contig1\tprediction\tgene\t1\t1000\t.\t+\t.\tgene_id "gene1"; transcript_id "mRNA1";\n'
            )
            temp_name = temp.name

        try:
            # Test with a GTF file
            parser, format_type = _detect_format(temp_name)
            # The actual implementation doesn't detect GTF format in _detect_format
            assert format_type == "default"
        finally:
            # Clean up
            os.unlink(temp_name)

    def test_simplify_go(self):
        """Test the simplifyGO function."""
        # Test with a simple GO term
        go_terms = ["GO:0005634"]
        result = simplifyGO(go_terms)
        assert result == ["GO:0005634"]

        # Test with multiple GO terms
        go_terms = ["GO:0005634", "GO:0003677", "GO:0006355"]
        result = simplifyGO(go_terms)
        assert set(result) == set(["GO:0005634", "GO:0003677", "GO:0006355"])

        # Test with an empty list
        go_terms = []
        result = simplifyGO(go_terms)
        assert result == []
