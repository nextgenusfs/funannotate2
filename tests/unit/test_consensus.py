"""
Unit tests for the consensus module.
"""

from gfftk.consensus import (
    auto_score_threshold,
    cluster_interlap,
    contained,
    ensure_unique_names,
    fasta_length,
    get_overlap,
)
from gfftk import interlap


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

    def test_cluster_interlap_no_duplicate_genes(self):
        """Test that cluster_interlap doesn't assign the same gene to multiple loci."""
        # Create an interlap object with overlapping gene models
        inter = interlap.InterLap()

        # Add gene models that reproduce the bug scenario
        # Gene models that should be in locus 1 (304791-308038)
        inter.add([304968, 307781, "gene1", "source1", [], 1])
        inter.add([307838, 307861, "gene2", "source2", [], 1])
        inter.add([308015, 308038, "gene3", "source3", [], 1])

        # Gene model that spans across loci (308249-309835) - this was causing the bug
        inter.add([308249, 309835, "problematic_gene", "helixer", [], 1])

        # Gene models that should be in locus 2 (308038-310019)
        inter.add([308954, 309835, "gene4", "source4", [], 1])
        inter.add([309000, 310019, "gene5", "source5", [], 1])

        # Cluster the genes
        clusters = cluster_interlap(inter)

        # Collect all gene IDs from all clusters
        all_gene_ids = []
        for cluster in clusters:
            for gene in cluster["genes"]:
                gene_id = gene[0]  # Gene ID is first element
                all_gene_ids.append(gene_id)

        # Check that no gene appears in multiple clusters
        unique_gene_ids = set(all_gene_ids)
        assert len(all_gene_ids) == len(unique_gene_ids), (
            f"Duplicate genes found: {all_gene_ids}"
        )

        # Verify that the problematic gene appears only once
        problematic_gene_count = all_gene_ids.count("problematic_gene")
        assert problematic_gene_count == 1, (
            f"Problematic gene appears {problematic_gene_count} times, should be 1"
        )

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
