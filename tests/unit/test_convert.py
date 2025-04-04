"""
Unit tests for the convert module.
"""
import os
import pytest
from gfftk.convert import (
    gff2bed,
    gff2gtf,
    gff2tbl,
)


class TestConvert:
    """Tests for conversion functions."""

    def test_gff2bed(self, sample_gff_file, temp_dir):
        """Test converting GFF to BED format."""
        # Set up output file
        output_file = os.path.join(temp_dir, "output.bed")

        # Convert GFF to BED
        gff2bed(sample_gff_file, output_file)

        # Check that the output file exists
        assert os.path.exists(output_file)

        # Check the content of the output file
        with open(output_file, "r") as f:
            content = f.readlines()

        # Basic checks on the content
        assert len(content) > 0

        # Check the format of the first line (should be BED format)
        first_line = content[0].strip().split("\t")
        assert len(first_line) >= 3  # BED format has at least 3 columns
        assert first_line[0] == "contig1"  # First column is chromosome/contig
        assert first_line[1].isdigit()  # Second column is start position
        assert first_line[2].isdigit()  # Third column is end position

    def test_gff2gtf(self, sample_gff_file, temp_dir):
        """Test converting GFF to GTF format."""
        # Set up output file
        output_file = os.path.join(temp_dir, "output.gtf")

        # Convert GFF to GTF
        gff2gtf(sample_gff_file, output_file)

        # Check that the output file exists
        assert os.path.exists(output_file)

        # Check the content of the output file
        with open(output_file, "r") as f:
            content = f.readlines()

        # Basic checks on the content
        assert len(content) > 0

        # Check for GTF-specific formatting in the attributes column
        for line in content:
            if line.startswith("#"):
                continue

            columns = line.strip().split("\t")
            if len(columns) < 9:
                continue

            attributes = columns[8]
            # GTF attributes are in the format 'key "value";'
            assert '"' in attributes
            assert ";" in attributes

    def test_gff2tbl(self, sample_gff_file, sample_fasta_file, temp_dir):
        """Test converting GFF to TBL format."""
        # Set up output file
        output_file = os.path.join(temp_dir, "output.tbl")

        # Convert GFF to TBL
        gff2tbl(sample_gff_file, sample_fasta_file, output_file)

        # Check that the output file exists
        assert os.path.exists(output_file)

        # Check the content of the output file
        with open(output_file, "r") as f:
            content = f.readlines()

        # Basic checks on the content
        assert len(content) > 0

        # Check for TBL-specific formatting
        # TBL format typically has '>Feature contig1' as a header
        assert any(line.startswith(">Feature") for line in content)
