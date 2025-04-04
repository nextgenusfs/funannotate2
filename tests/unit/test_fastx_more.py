"""
More unit tests for the fastx module.
"""

import os
import pytest
from funannotate2.fastx import (
    contig_analysis,
    analyzeAssembly,
    analyzeAssemblySimple,
)


class TestContigAnalysis:
    """Tests for the contig_analysis function."""

    def test_contig_analysis(self):
        """Test analyzing a sequence for masked regions and gaps."""
        # Create a test sequence with masked regions and gaps
        title = "test_sequence"
        seq = "ATGCGTacgtacgtNNNNNatgcATGC"

        # Run the function
        result_title, masked, gaps = contig_analysis(title, seq)

        # Check the results
        assert result_title == title

        # Check masked regions (lowercase and N's)
        # The function groups consecutive masked regions together
        assert len(masked) == 1
        assert (6, 22) in masked  # combined lowercase and N's region

        # Check gaps (only N's)
        assert len(gaps) == 1
        assert (14, 18) in gaps  # 0-based indexing


class TestAnalyzeAssembly:
    """Tests for the analyzeAssembly function."""

    def test_analyze_assembly(self, temp_dir, test_data_dir):
        """Test analyzing a valid FASTA file."""
        # Setup input and output files
        fasta_file = os.path.join(test_data_dir, "test_valid.fasta")
        masked_output = os.path.join(temp_dir, "masked.bed")
        gaps_output = os.path.join(temp_dir, "gaps.bed")

        # Run the function
        stats, bad_names, errors, _ = analyzeAssembly(
            fasta_file, masked_output, gaps_output, header_max=50
        )

        # Check the results
        assert isinstance(stats, dict)
        assert "n_contigs" in stats
        assert stats["n_contigs"] == 3
        assert "size" in stats

        # No bad names or errors expected
        assert bad_names == []
        assert errors == []

        # Check that output files were created
        assert os.path.exists(masked_output)
        assert os.path.exists(gaps_output)


class TestAnalyzeAssemblySimple:
    """Tests for the analyzeAssemblySimple function."""

    def test_analyze_assembly_simple_valid(self, test_data_dir):
        """Test getting stats for a valid FASTA file."""
        # Setup input file
        fasta_file = os.path.join(test_data_dir, "test_valid.fasta")

        # Run the function
        stats, bad_names, errors = analyzeAssemblySimple(fasta_file)

        # Check the results
        assert isinstance(stats, dict)
        assert "n_contigs" in stats
        assert stats["n_contigs"] == 3
        assert "size" in stats
        assert "n50" in stats
        assert "n90" in stats
        assert "l50" in stats
        assert "l90" in stats
        assert "avg_length" in stats

        # No bad names or errors expected
        assert bad_names == []
        assert errors == []

    def test_analyze_assembly_simple_invalid(self, test_data_dir):
        """Test getting stats for an invalid FASTA file."""
        # Setup input file with invalid characters and long headers
        fasta_file = os.path.join(test_data_dir, "test_invalid.fasta")

        # Run the function
        stats, bad_names, errors = analyzeAssemblySimple(fasta_file, header_max=50)

        # Check the results
        assert isinstance(stats, dict)
        assert "n_contigs" in stats
        assert stats["n_contigs"] == 2

        # Check for bad names (long headers)
        assert len(bad_names) == 1
        assert "contig2_with_very_long_header" in bad_names[0]

        # Check for errors (invalid characters)
        assert len(errors) == 1
        assert errors[0][0] == "X"  # The invalid character is X
