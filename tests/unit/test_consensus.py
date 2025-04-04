"""
Unit tests for the consensus module.
"""
import os
import pytest
from gfftk.consensus import (
    get_overlap,
    contained,
    auto_score_threshold,
    ensure_unique_names,
    fasta_length,
)


class TestConsensusHelpers:
    """Tests for helper functions in the consensus module."""

    def test_get_overlap(self):
        """Test the get_overlap function."""
        # Test with overlapping ranges
        a = [10, 20]
        b = [15, 25]
        overlap = get_overlap(a, b)
        assert overlap == 5  # Overlap is 15-20 = 5

        # Test with non-overlapping ranges
        a = [10, 20]
        b = [30, 40]
        overlap = get_overlap(a, b)
        assert overlap == 0

        # Test with one range contained in the other
        a = [10, 30]
        b = [15, 25]
        overlap = get_overlap(a, b)
        assert overlap == 10  # Overlap is 15-25 = 10

        # Test with ranges touching but not overlapping
        a = [10, 20]
        b = [20, 30]
        overlap = get_overlap(a, b)
        assert overlap == 0

    def test_contained(self):
        """Test the contained function."""
        # Test with one range contained in the other
        a = [15, 25]
        b = [10, 30]
        result = contained(a, b)
        assert result is True

        # Test with ranges that overlap but neither contains the other
        a = [10, 20]
        b = [15, 25]
        result = contained(a, b)
        assert result is False

        # Test with non-overlapping ranges
        a = [10, 20]
        b = [30, 40]
        result = contained(a, b)
        assert result is False

        # Test with identical ranges
        a = [10, 20]
        b = [10, 20]
        result = contained(a, b)
        assert result is True

    def test_auto_score_threshold(self):
        """Test the auto_score_threshold function."""
        # Test with default weights
        weights = {"source1": 1, "source2": 2, "source3": 3}
        order = ["source1", "source2", "source3"]
        threshold = auto_score_threshold(weights, order)
        assert threshold == 6  # Default user_weight is 6

        # Test with custom user_weight
        threshold = auto_score_threshold(weights, order, user_weight=10)
        assert threshold == 10

        # Test with empty weights
        weights = {}
        order = []
        threshold = auto_score_threshold(weights, order)
        assert threshold == 6  # Default user_weight is 6

    def test_ensure_unique_names(self):
        """Test the ensure_unique_names function."""
        # Create a sample gene dictionary with duplicate names
        genes = {
            "gene1": {"name": "duplicate_name"},
            "gene2": {"name": "duplicate_name"},
            "gene3": {"name": "unique_name"},
        }

        # Ensure unique names
        result = ensure_unique_names(genes)

        # Check that the names are now unique
        names = [gene["name"] for gene in result.values()]
        assert len(names) == len(set(names))

        # Check that the original unique name is preserved
        assert result["gene3"]["name"] == "unique_name"

        # Check that the duplicate names are modified
        assert result["gene1"]["name"] != result["gene2"]["name"]
        assert result["gene1"]["name"].startswith("duplicate_name")
        assert result["gene2"]["name"].startswith("duplicate_name")

    def test_fasta_length(self, sample_fasta_file):
        """Test the fasta_length function."""
        # Get the length of the FASTA file
        length = fasta_length(sample_fasta_file)

        # Check that the length is correct
        assert length == 180  # 3 lines of 60 characters each
