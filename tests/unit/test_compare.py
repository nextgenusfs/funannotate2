"""
Unit tests for the compare module.
"""
import os
import tempfile
from gfftk.compare import (
    compare_gff,
    calculate_stats,
    calculate_sensitivity_specificity,
)


class TestCompare:
    """Tests for compare functions."""

    def test_calculate_sensitivity_specificity(self):
        """Test calculating sensitivity and specificity."""
        # Test case 1: Perfect prediction
        tp = 100  # True positives
        fp = 0  # False positives
        fn = 0  # False negatives

        result = calculate_sensitivity_specificity(tp, fp, fn)

        assert result["sensitivity"] == 100.0
        assert result["specificity"] == 100.0
        assert result["precision"] == 100.0
        assert result["f1"] == 100.0

        # Test case 2: Imperfect prediction
        tp = 80  # True positives
        fp = 20  # False positives
        fn = 20  # False negatives

        result = calculate_sensitivity_specificity(tp, fp, fn)

        assert result["sensitivity"] == 80.0
        assert result["specificity"] == 80.0
        assert result["precision"] == 80.0
        assert result["f1"] == 80.0

        # Test case 3: No true positives
        tp = 0  # True positives
        fp = 50  # False positives
        fn = 50  # False negatives

        result = calculate_sensitivity_specificity(tp, fp, fn)

        assert result["sensitivity"] == 0.0
        assert result["specificity"] == 0.0
        assert result["precision"] == 0.0
        assert result["f1"] == 0.0

    def test_calculate_stats(self):
        """Test calculating comparison statistics."""
        # Create sample data
        ref_genes = {
            "gene1": {
                "type": "gene",
                "contig": "contig1",
                "location": [1, 1000],
                "strand": "+",
            },
            "gene2": {
                "type": "gene",
                "contig": "contig1",
                "location": [2000, 3000],
                "strand": "-",
            },
            "gene3": {
                "type": "gene",
                "contig": "contig2",
                "location": [1, 1000],
                "strand": "+",
            },
        }

        pred_genes = {
            "pred1": {
                "type": "gene",
                "contig": "contig1",
                "location": [1, 1000],
                "strand": "+",
            },
            "pred2": {
                "type": "gene",
                "contig": "contig1",
                "location": [2100, 3100],
                "strand": "-",
            },
            "pred4": {
                "type": "gene",
                "contig": "contig2",
                "location": [2000, 3000],
                "strand": "+",
            },
        }

        # Calculate stats
        result = calculate_stats(ref_genes, pred_genes)

        # Check the result
        assert "gene" in result
        assert result["gene"]["total_ref"] == 3
        assert result["gene"]["total_pred"] == 3
        assert result["gene"]["found"] == 1  # Only gene1/pred1 is an exact match

        # Check sensitivity/specificity
        assert "sensitivity" in result["gene"]
        assert "specificity" in result["gene"]
        assert "precision" in result["gene"]
        assert "f1" in result["gene"]

    def test_compare_gff_basic(self):
        """Test basic GFF comparison."""
        # Create a reference GFF file
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".gff3"
        ) as ref_gff:
            ref_gff.write("##gff-version 3\n")
            ref_gff.write(
                "contig1\tprediction\tgene\t1\t1000\t.\t+\t.\tID=gene1;Name=test_gene1\n"
            )
            ref_gff.write(
                "contig1\tprediction\tgene\t2000\t3000\t.\t-\t.\tID=gene2;Name=test_gene2\n"
            )
            ref_gff_name = ref_gff.name

        # Create a prediction GFF file
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".gff3"
        ) as pred_gff:
            pred_gff.write("##gff-version 3\n")
            pred_gff.write(
                "contig1\tprediction\tgene\t1\t1000\t.\t+\t.\tID=pred1;Name=pred_gene1\n"
            )
            pred_gff.write(
                "contig1\tprediction\tgene\t4000\t5000\t.\t+\t.\tID=pred2;Name=pred_gene2\n"
            )
            pred_gff_name = pred_gff.name

        # Create a temporary FASTA file
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".fasta"
        ) as temp_fasta:
            temp_fasta.write(">contig1\n")
            temp_fasta.write("A" * 5000 + "\n")
            temp_fasta_name = temp_fasta.name

        try:
            # Compare the GFF files
            result = compare_gff(ref_gff_name, pred_gff_name, temp_fasta_name)

            # Check the result
            assert "stats" in result
            assert "gene" in result["stats"]
            assert result["stats"]["gene"]["total_ref"] == 2
            assert result["stats"]["gene"]["total_pred"] == 2
            assert (
                result["stats"]["gene"]["found"] == 1
            )  # Only gene1/pred1 is an exact match
        finally:
            # Clean up
            for filename in [ref_gff_name, pred_gff_name, temp_fasta_name]:
                if os.path.exists(filename):
                    os.unlink(filename)
