"""
Unit tests for the paf module.
"""
import os
import tempfile
import pytest
from gfftk.paf import (
    paf2dict,
    cs2coords,
    cs2tuples,
)


class TestPAF:
    """Tests for PAF functions."""

    def test_cs2tuples(self):
        """Test parsing CIGAR strings into tuples."""
        # Test with a simple CIGAR string
        cs = ":60~gt63ag:365"
        result = cs2tuples(cs)

        # Check the result
        assert len(result) == 3
        assert result[0] == (":", "60")
        assert result[1] == ("~", "gt63ag")
        assert result[2] == (":", "365")

        # Test with a more complex CIGAR string
        cs = ":60~gt63ag:365~gt49ag:520-atgc+cta:100"
        result = cs2tuples(cs)

        # Check the result
        assert len(result) == 6
        assert result[0] == (":", "60")
        assert result[1] == ("~", "gt63ag")
        assert result[2] == (":", "365")
        assert result[3] == ("~", "gt49ag")
        assert result[4] == (":", "520")
        assert result[5] == ("-", "atgc")

        # Test with an empty string
        cs = ""
        result = cs2tuples(cs)
        assert len(result) == 0

    def test_cs2coords_plus_strand(self):
        """Test converting CIGAR string to coordinates on plus strand."""
        # Test with a simple CIGAR string on plus strand
        start = 1000
        qstart = 0
        length = 100
        strand = "+"
        cs = "cs:Z::60~gt63ag:40"

        result = cs2coords(start, qstart, length, strand, cs)

        # Check the result
        exons, query, proper_splice = result

        # Should have 2 exons
        assert len(exons) == 2

        # First exon should start at 1001 (start + offset) and end at 1060
        assert exons[0] == (1001, 1060)

        # Second exon should start at 1124 (1060 + 63 + 1) and end at 1164 (1124 + 40)
        assert exons[1] == (1124, 1164)

        # Query coordinates should match
        assert len(query) == 2
        assert query[0] == (0, 60)
        assert query[1] == (60, 100)

        # Should have proper splice sites
        assert proper_splice is True

    def test_cs2coords_minus_strand(self):
        """Test converting CIGAR string to coordinates on minus strand."""
        # Test with a simple CIGAR string on minus strand
        start = 1000
        qstart = 0
        length = 100
        strand = "-"
        cs = "cs:Z::60~ct63ac:40"

        result = cs2coords(start, qstart, length, strand, cs)

        # Check the result
        exons, query, proper_splice = result

        # Should have 2 exons
        assert len(exons) == 2

        # Exons should be in reverse order for minus strand
        # First exon should start at 1124 and end at 1164
        # Second exon should start at 1001 and end at 1060
        assert exons[0] == (1124, 1164)
        assert exons[1] == (1001, 1060)

        # Query coordinates should match
        assert len(query) == 2
        assert query[0] == (0, 60)
        assert query[1] == (60, 100)

        # Should have proper splice sites
        assert proper_splice is True

    def test_cs2coords_with_indels(self):
        """Test converting CIGAR string with indels to coordinates."""
        # Test with a CIGAR string containing indels
        start = 1000
        qstart = 0
        length = 100
        strand = "+"
        cs = "cs:Z::50-atgc:10+cta:40"

        result = cs2coords(start, qstart, length, strand, cs)

        # Check the result
        exons, query, proper_splice = result

        # Should have 1 exon
        assert len(exons) == 1

        # Exon should start at 1001 and end at 1104 (1001 + 50 + 4 + 10 + 3 + 40 - 4)
        # The -4 is because we have a deletion of 4 bases
        assert exons[0] == (1001, 1104)

        # Query coordinates should match
        assert len(query) == 1
        assert query[0] == (0, 100)

        # No splice sites, so proper_splice should be True
        assert proper_splice is True

    def test_cs2coords_improper_splice(self):
        """Test converting CIGAR string with improper splice sites."""
        # Test with a CIGAR string containing improper splice sites
        start = 1000
        qstart = 0
        length = 100
        strand = "+"
        cs = "cs:Z::60~xx63yy:40"  # Using xx and yy instead of gt and ag

        result = cs2coords(start, qstart, length, strand, cs)

        # Check the result
        exons, query, proper_splice = result

        # Should have 2 exons
        assert len(exons) == 2

        # First exon should start at 1001 and end at 1060
        assert exons[0] == (1001, 1060)

        # Second exon should start at 1124 and end at 1164
        assert exons[1] == (1124, 1164)

        # Query coordinates should match
        assert len(query) == 2
        assert query[0] == (0, 60)
        assert query[1] == (60, 100)

        # Should have improper splice sites
        assert proper_splice is False

    def test_paf2dict(self):
        """Test parsing PAF file to dictionary."""
        # Get the path to the test data
        test_dir = os.path.dirname(os.path.abspath(__file__))
        paf_file = os.path.join(test_dir, "test_data", "sample.paf")
        fasta_file = os.path.join(test_dir, "test_data", "sample.fasta")

        # Parse the PAF file
        result = paf2dict(paf_file, fasta_file)

        # Check the result
        assert isinstance(result, dict)
        assert len(result) == 2

        # Check the first entry
        assert "OPO1_006208-T1" in result
        entry = result["OPO1_006208-T1"]
        assert len(entry) == 4

        # Check the exons
        exons, query, proper_splice, cs = entry
        assert len(exons) == 9  # 9 exons in the first entry
        assert proper_splice is True
        assert cs.startswith("cs:Z:")

        # Check the second entry
        assert "OPO1_006778-T1" in result
        entry = result["OPO1_006778-T1"]
        assert len(entry) == 4

        # Check the exons
        exons, query, proper_splice, cs = entry
        assert len(exons) == 8  # 8 exons in the second entry
        assert proper_splice is True
        assert cs.startswith("cs:Z:")

    def test_paf2dict_with_min_mapq(self):
        """Test parsing PAF file with minimum mapping quality."""
        # Create a temporary PAF file with different mapping qualities
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".paf"
        ) as temp_paf:
            # Entry with high mapping quality
            temp_paf.write(
                "OPO1_006208-T1\t2883\t0\t2883\t+\tscaffold_60\t41779\t19205\t22545\t2883\t2883\t60\tcs:Z::60~gt63ag:365\n"
            )
            # Entry with low mapping quality
            temp_paf.write(
                "OPO1_006778-T1\t3150\t0\t3150\t-\tscaffold_82\t14472\t8938\t12444\t3150\t3152\t1\tcs:Z::377~ct45ac:171\n"
            )
            temp_paf_name = temp_paf.name

        # Create a temporary FASTA file
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".fasta"
        ) as temp_fasta:
            temp_fasta.write(">scaffold_60\nATGCATGCATGC\n")
            temp_fasta.write(">scaffold_82\nGCTAGCTAGCTA\n")
            temp_fasta_name = temp_fasta.name

        try:
            # Parse the PAF file with default min_mapq (2)
            result = paf2dict(temp_paf_name, temp_fasta_name)

            # Check the result
            assert isinstance(result, dict)
            assert len(result) == 1  # Only the high-quality entry should be included
            assert "OPO1_006208-T1" in result
            assert "OPO1_006778-T1" not in result

            # Parse the PAF file with min_mapq=0
            result = paf2dict(temp_paf_name, temp_fasta_name, min_mapq=0)

            # Check the result
            assert isinstance(result, dict)
            assert len(result) == 2  # Both entries should be included
            assert "OPO1_006208-T1" in result
            assert "OPO1_006778-T1" in result
        finally:
            # Clean up
            os.unlink(temp_paf_name)
            os.unlink(temp_fasta_name)

    def test_paf2dict_with_annotation(self):
        """Test parsing PAF file with existing annotation."""
        # Get the path to the test data
        test_dir = os.path.dirname(os.path.abspath(__file__))
        paf_file = os.path.join(test_dir, "test_data", "sample.paf")
        fasta_file = os.path.join(test_dir, "test_data", "sample.fasta")

        # Create an existing annotation
        annotation = {"existing_entry": ("exons", "query", True, "cs:Z::100")}

        # Parse the PAF file with the existing annotation
        result = paf2dict(paf_file, fasta_file, annotation=annotation)

        # Check the result
        assert isinstance(result, dict)
        assert len(result) == 3  # 2 from the PAF file + 1 from the existing annotation

        # Check that the existing entry is preserved
        assert "existing_entry" in result
        assert result["existing_entry"] == ("exons", "query", True, "cs:Z::100")

        # Check that the new entries are added
        assert "OPO1_006208-T1" in result
        assert "OPO1_006778-T1" in result

    def test_paf2dict_without_cs_tag(self):
        """Test parsing PAF file without cs tag."""
        # Create a temporary PAF file without cs tag
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".paf"
        ) as temp_paf:
            temp_paf.write(
                "OPO1_006208-T1\t2883\t0\t2883\t+\tscaffold_60\t41779\t19205\t22545\t2883\t2883\t60\n"
            )
            temp_paf_name = temp_paf.name

        # Create a temporary FASTA file
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".fasta"
        ) as temp_fasta:
            temp_fasta.write(">scaffold_60\nATGCATGCATGC\n")
            temp_fasta_name = temp_fasta.name

        try:
            # Parse the PAF file
            result = paf2dict(temp_paf_name, temp_fasta_name)

            # Check the result
            assert isinstance(result, dict)
            assert (
                len(result) == 0
            )  # No entries should be included because there's no cs tag
        finally:
            # Clean up
            os.unlink(temp_paf_name)
            os.unlink(temp_fasta_name)
