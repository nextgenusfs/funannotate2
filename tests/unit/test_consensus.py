"""
Unit tests for the consensus module.
"""

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
        # The actual implementation returns False for identical ranges
        # because left = 0 and right = 0, so left > 0 and right < 0 is False
        a = [10, 20]
        b = [10, 20]
        result = contained(a, b)
        assert result is False

    def test_auto_score_threshold(self):
        """Test the auto_score_threshold function."""
        # Test with default weights
        weights = {"source1": 1, "source2": 2, "source3": 3}
        # The actual implementation expects order to be a dictionary, not a list
        order = {"source1": 1, "source2": 2, "source3": 3}
        threshold = auto_score_threshold(weights, order)
        # The actual implementation calculates threshold as 1 + min(allweights.values())
        # where allweights[w] = weights.get(w, 1) * user_weight
        # So the minimum weight is 1 * 6 = 6, and threshold is 1 + 6 = 7
        assert threshold == 7

        # Test with custom user_weight
        threshold = auto_score_threshold(weights, order, user_weight=10)
        # The minimum weight is 1 * 10 = 10, and threshold is 1 + 10 = 11
        assert threshold == 11

        # Test with empty weights but non-empty order
        # The actual implementation will raise a ValueError if order is empty
        # because it tries to calculate min() on an empty sequence
        weights = {}
        order = {"source1": 1}  # Need at least one item in order
        threshold = auto_score_threshold(weights, order)
        # With empty weights, the function will use the default weight of 1
        # So the minimum weight is 1 * 6 = 6, and threshold is 1 + 6 = 7
        assert threshold == 7

    def test_ensure_unique_names(self):
        """Test the ensure_unique_names function."""
        # Create a sample gene dictionary
        genes = {
            "gene1": {"name": "gene1_name"},
            "gene2": {"name": "gene2_name"},
            "gene3": {"name": "gene3_name"},
        }

        # Ensure unique names
        result = ensure_unique_names(genes)

        # The actual implementation doesn't modify the gene names
        # It adds a unique slug to the gene IDs
        assert len(result) == 3

        # Check that the original gene data is preserved
        for key in result:
            # The key should be in the format "gene{n}.{slug}"
            assert key.startswith("gene")
            assert "." in key

            # The gene data should be the same as the original
            gene_id = key.split(".")[0]
            assert result[key] == genes[gene_id]

    def test_fasta_length(self, sample_fasta_file):
        """Test the fasta_length function."""
        # Get the length of the FASTA file
        length = fasta_length(sample_fasta_file)

        # Check that the length is correct
        assert length == 180  # 3 lines of 60 characters each
