"""
Unit tests for the fasta module.
"""
import os
import tempfile
from gfftk.fasta import (
    fasta2dict,
    dict2fasta,
    fasta2headers,
    fasta_stats,
    RevComp,
    translate,
    getSeqRegions,
)


class TestFasta:
    """Tests for FASTA functions."""

    def test_fasta2dict(self):
        """Test converting FASTA to dictionary."""
        # Create a temporary FASTA file
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".fasta"
        ) as temp:
            temp.write(">seq1\n")
            temp.write("ATGCATGCATGC\n")
            temp.write(">seq2\n")
            temp.write("GTATCGATCGAT\n")
            temp_name = temp.name

        try:
            # Convert FASTA to dictionary
            result = fasta2dict(temp_name)

            # Check the result
            assert len(result) == 2
            assert "seq1" in result
            assert "seq2" in result
            assert result["seq1"] == "ATGCATGCATGC"
            assert result["seq2"] == "GTATCGATCGAT"
        finally:
            # Clean up
            os.unlink(temp_name)

    def test_dict2fasta(self):
        """Test converting dictionary to FASTA."""
        # Create a dictionary
        fasta_dict = {
            "seq1": "ATGCATGCATGC",
            "seq2": "GTATCGATCGAT",
        }

        # Create a temporary output file
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".fasta"
        ) as temp:
            temp_name = temp.name

        try:
            # Convert dictionary to FASTA
            dict2fasta(fasta_dict, temp_name)

            # Check that the output file exists
            assert os.path.exists(temp_name)

            # Check the content of the output file
            with open(temp_name, "r") as f:
                content = f.read()

            # Basic checks on the content
            assert ">seq1" in content
            assert "ATGCATGCATGC" in content
            assert ">seq2" in content
            assert "GTATCGATCGAT" in content
        finally:
            # Clean up
            os.unlink(temp_name)

    def test_fasta2headers(self):
        """Test extracting headers from FASTA."""
        # Create a temporary FASTA file
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".fasta"
        ) as temp:
            temp.write(">seq1 description 1\n")
            temp.write("ATGCATGCATGC\n")
            temp.write(">seq2 description 2\n")
            temp.write("GTATCGATCGAT\n")
            temp_name = temp.name

        try:
            # Extract headers from FASTA
            result = fasta2headers(temp_name)

            # Check the result
            assert len(result) == 2
            assert "seq1" in result
            assert "seq2" in result
            assert result["seq1"] == "description 1"
            assert result["seq2"] == "description 2"
        finally:
            # Clean up
            os.unlink(temp_name)

    def test_fasta_stats(self):
        """Test calculating FASTA statistics."""
        # Create a temporary FASTA file
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".fasta"
        ) as temp:
            temp.write(">seq1\n")
            temp.write("ATGCATGCATGC\n")  # 12 bp
            temp.write(">seq2\n")
            temp.write("GTATCGATCGAT\n")  # 12 bp
            temp_name = temp.name

        try:
            # Calculate FASTA statistics
            result = fasta_stats(temp_name)

            # Check the result
            assert result["records"] == 2
            assert result["sum"] == 24
            assert result["mean"] == 12.0
            assert result["min"] == 12
            assert result["max"] == 12
            assert result["N50"] == 12
            assert result["GC"] == 0.5  # 12 G/C out of 24 total bases
        finally:
            # Clean up
            os.unlink(temp_name)

    def test_RevComp(self):
        """Test reverse complementing DNA sequences."""
        # Test with a simple sequence
        seq = "ATGCATGCATGC"
        result = RevComp(seq)
        assert result == "GCATGCATGCAT"

        # Test with a sequence containing N
        seq = "ATGCNNGCATGC"
        result = RevComp(seq)
        assert result == "GCATGCNNGCAT"

        # Test with an empty sequence
        seq = ""
        result = RevComp(seq)
        assert result == ""

    def test_translate(self):
        """Test translating DNA to protein."""
        # Test with a simple sequence
        seq = "ATGCATGCATGC"
        result = translate(seq)
        assert result == "MHAC"

        # Test with a sequence containing N
        seq = "ATGCNNGCATGC"
        result = translate(seq)
        assert "X" in result  # N translates to X

        # Test with a sequence that's not a multiple of 3
        seq = "ATGCATGCATG"  # 11 bp
        result = translate(seq)
        assert len(result) == 3  # 11 // 3 = 3 complete codons

        # Test with an empty sequence
        seq = ""
        result = translate(seq)
        assert result == ""

    def test_getSeqRegions(self):
        """Test extracting regions from sequences."""
        # Create a dictionary of sequences
        fasta_dict = {
            "seq1": "ATGCATGCATGC",  # 12 bp
            "seq2": "GTATCGATCGAT",  # 12 bp
        }

        # Define regions to extract
        regions = [
            ["seq1", 1, 6],  # ATGCAT
            ["seq1", 7, 12],  # GCATGC
            ["seq2", 4, 9],  # CGATCG
        ]

        # Extract regions
        result = getSeqRegions(fasta_dict, regions)

        # Check the result
        assert len(result) == 3
        assert result[0] == "ATGCAT"
        assert result[1] == "GCATGC"
        assert result[2] == "CGATCG"

        # Test with regions that are out of bounds
        regions = [
            ["seq1", 1, 20],  # Out of bounds
            ["seq3", 1, 6],  # Non-existent sequence
        ]

        # Extract regions
        result = getSeqRegions(fasta_dict, regions)

        # Check the result
        assert len(result) == 2
        assert result[0] == ""  # Out of bounds returns empty string
        assert result[1] == ""  # Non-existent sequence returns empty string
