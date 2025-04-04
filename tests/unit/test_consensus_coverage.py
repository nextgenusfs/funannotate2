"""
Unit tests for the enhanced consensus module with coverage-based scoring.
"""

from gfftk.consensus import (
    score_evidence,
    map_coords,
)


class TestConsensusCoverage:
    """Tests for consensus functions with coverage-based scoring."""

    def test_score_evidence_exact_match(self):
        """Test the score_evidence function with exact match."""
        # Test with identical coordinates
        g_coords = [[1, 100], [200, 300]]
        e_coords = [[1, 100], [200, 300]]

        score = score_evidence(g_coords, e_coords)
        assert score == 20  # Perfect match with default weight of 2 (10 * 2)

    def test_score_evidence_partial_coverage(self):
        """Test the score_evidence function with partial coverage."""
        # Test with partial coverage (50%)
        g_coords = [[1, 100], [200, 300]]  # Total length: 200
        e_coords = [[25, 75], [200, 250]]  # Covers: 100 (50%)

        score = score_evidence(g_coords, e_coords)
        # Base score would be 7.5 (average of 5 and 10), adjusted by coverage (50%)
        # 7.5 * (0.5 + 0.5 * 0.5) = 7.5 * 0.75 = 5.625, rounded to 5 * 2 = 10
        assert score == 10

    def test_score_evidence_high_coverage(self):
        """Test the score_evidence function with high coverage."""
        # Test with high coverage (80%)
        g_coords = [[1, 100], [200, 300]]  # Total length: 200
        e_coords = [[10, 100], [200, 280]]  # Covers: 170 (85%)

        score = score_evidence(g_coords, e_coords)
        # Base score would be 7.5, adjusted by coverage (85%)
        # 7.5 * (0.5 + 0.5 * 0.85) = 7.5 * 0.925 = 6.9375, rounded to 7 * 2 = 14
        assert score == 14

    def test_score_evidence_low_coverage(self):
        """Test the score_evidence function with low coverage."""
        # Test with low coverage (20%)
        g_coords = [[1, 100], [200, 300]]  # Total length: 200
        e_coords = [[1, 20], [200, 220]]  # Covers: 40 (20%)

        score = score_evidence(g_coords, e_coords)
        # Base score would be 7.5, adjusted by coverage (20%)
        # 7.5 * (0.5 + 0.5 * 0.2) = 7.5 * 0.6 = 4.5, rounded to 4 * 2 = 8
        assert score == 8

    def test_score_evidence_no_overlap(self):
        """Test the score_evidence function with no overlap."""
        # Test with no overlap
        g_coords = [[1, 100], [200, 300]]
        e_coords = [[400, 500], [600, 700]]

        score = score_evidence(g_coords, e_coords)
        assert score == 0  # No overlap should have score of 0

    def test_score_evidence_single_exon(self):
        """Test the score_evidence function with a single exon gene."""
        # Test with a single exon gene, perfect match
        g_coords = [[1, 100]]
        e_coords = [[1, 100]]

        score = score_evidence(g_coords, e_coords)
        assert score == 20  # Perfect match with default weight of 2 (10 * 2)

        # Test with a single exon gene, partial coverage (50%)
        g_coords = [[1, 100]]
        e_coords = [[25, 75]]

        score = score_evidence(g_coords, e_coords)
        # Base score would be 10, adjusted by coverage (50%)
        # 10 * (0.5 + 0.5 * 0.5) = 10 * 0.75 = 7.5, rounded to 7 * 2 = 14
        assert score == 14

    def test_score_evidence_custom_weight(self):
        """Test the score_evidence function with custom weight."""
        # Test with custom weight
        g_coords = [[1, 100], [200, 300]]
        e_coords = [[1, 100], [200, 300]]

        score = score_evidence(g_coords, e_coords, weight=5)
        assert score == 50  # Perfect match with weight of 5 (10 * 5)

        # Test with partial coverage and custom weight
        g_coords = [[1, 100], [200, 300]]  # Total length: 200
        e_coords = [[25, 75], [200, 250]]  # Covers: 100 (50%)

        score = score_evidence(g_coords, e_coords, weight=5)
        # Base score would be 7.5, adjusted by coverage (50%)
        # 7.5 * (0.5 + 0.5 * 0.5) = 7.5 * 0.75 = 5.625, rounded to 5 * 5 = 25
        assert score == 25
