"""
Unit tests for the filter_and_write_gff3 function in utilities.py.
"""

import os
import tempfile
import pytest

from funannotate2.utilities import filter_and_write_gff3


class TestFilterAndWriteGff3:
    """Tests for the filter_and_write_gff3 function."""

    @pytest.fixture
    def temp_genome_fasta(self):
        """Create a temporary genome FASTA file."""
        f = tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False)
        f.write(">contig1\n")
        f.write("ATGAAACGTAAGGCCTAG" * 1000)  # ~18000 bp
        f.write("\n")
        f.close()
        yield f.name
        os.unlink(f.name)

    @pytest.fixture
    def temp_gff3_file(self):
        """Create a temporary GFF3 file with test data."""
        f = tempfile.NamedTemporaryFile(mode="w", suffix=".gff3", delete=False)
        # Gene 1: 300 bp CDS = 100 aa (should pass with default min=30, max=30000)
        f.write("contig1\ttest\tgene\t1\t300\t.\t+\t.\tID=gene1\n")
        f.write("contig1\ttest\tmRNA\t1\t300\t.\t+\t.\tID=mRNA1;Parent=gene1\n")
        f.write("contig1\ttest\tCDS\t1\t300\t.\t+\t0\tID=CDS1;Parent=mRNA1\n")
        f.write("contig1\ttest\texon\t1\t300\t.\t+\t.\tID=exon1;Parent=mRNA1\n")

        # Gene 2: 90 bp CDS = 30 aa (should pass with default min=30)
        f.write("contig1\ttest\tgene\t400\t500\t.\t+\t.\tID=gene2\n")
        f.write("contig1\ttest\tmRNA\t400\t500\t.\t+\t.\tID=mRNA2;Parent=gene2\n")
        f.write("contig1\ttest\tCDS\t400\t500\t.\t+\t0\tID=CDS2;Parent=mRNA2\n")
        f.write("contig1\ttest\texon\t400\t500\t.\t+\t.\tID=exon2;Parent=mRNA2\n")

        # Gene 3: 60 bp CDS = 20 aa (should be filtered with default min=30)
        f.write("contig1\ttest\tgene\t600\t660\t.\t+\t.\tID=gene3\n")
        f.write("contig1\ttest\tmRNA\t600\t660\t.\t+\t.\tID=mRNA3;Parent=gene3\n")
        f.write("contig1\ttest\tCDS\t600\t660\t.\t+\t0\tID=CDS3;Parent=mRNA3\n")
        f.write("contig1\ttest\texon\t600\t660\t.\t+\t.\tID=exon3;Parent=mRNA3\n")

        # Gene 4: 90000 bp CDS = 30000 aa (should pass with default max=30000)
        f.write("contig1\ttest\tgene\t1000\t91000\t.\t+\t.\tID=gene4\n")
        f.write("contig1\ttest\tmRNA\t1000\t91000\t.\t+\t.\tID=mRNA4;Parent=gene4\n")
        f.write("contig1\ttest\tCDS\t1000\t91000\t.\t+\t0\tID=CDS4;Parent=mRNA4\n")
        f.write("contig1\ttest\texon\t1000\t91000\t.\t+\t.\tID=exon4;Parent=mRNA4\n")

        # Gene 5: 90003 bp CDS = 30001 aa (should be filtered with default max=30000)
        f.write("contig1\ttest\tgene\t100000\t190003\t.\t+\t.\tID=gene5\n")
        f.write("contig1\ttest\tmRNA\t100000\t190003\t.\t+\t.\tID=mRNA5;Parent=gene5\n")
        f.write("contig1\ttest\tCDS\t100000\t190003\t.\t+\t0\tID=CDS5;Parent=mRNA5\n")
        f.write("contig1\ttest\texon\t100000\t190003\t.\t+\t.\tID=exon5;Parent=mRNA5\n")
        f.close()
        yield f.name
        os.unlink(f.name)

    def test_filter_default_parameters(self, temp_gff3_file, temp_genome_fasta):
        """Test filtering with default parameters (min=30, max=30000)."""
        output_file = tempfile.mktemp(suffix=".gff3")

        try:
            stats = filter_and_write_gff3(temp_gff3_file, output_file)

            # Should keep genes 1, 2, 4 (3 genes)
            # Should filter genes 3 (too small), 5 (too large)
            assert stats["kept"] == 3
            assert stats["filtered"] == 2

            # Verify output file contains expected genes
            with open(output_file, "r") as f:
                content = f.read()
                assert "gene1" in content
                assert "gene2" in content
                assert "gene3" not in content
                assert "gene4" in content
                assert "gene5" not in content
        finally:
            if os.path.exists(output_file):
                os.unlink(output_file)

    def test_filter_custom_min_length(self, temp_gff3_file, temp_genome_fasta):
        """Test filtering with custom minimum protein length."""
        output_file = tempfile.mktemp(suffix=".gff3")

        try:
            # Set min to 50 aa (150 bp CDS)
            stats = filter_and_write_gff3(
                temp_gff3_file,
                output_file,
                min_protein_length=50,
                max_protein_length=30000,
            )

            # Should keep genes 1, 4 (2 genes)
            # Should filter genes 2 (too small), 3 (too small), 5 (too large)
            assert stats["kept"] == 2
            assert stats["filtered"] == 3

            # Verify output file
            with open(output_file, "r") as f:
                content = f.read()
                assert "gene1" in content
                assert "gene2" not in content
                assert "gene4" in content
        finally:
            if os.path.exists(output_file):
                os.unlink(output_file)

    def test_filter_custom_max_length(self, temp_gff3_file, temp_genome_fasta):
        """Test filtering with custom maximum protein length."""
        output_file = tempfile.mktemp(suffix=".gff3")

        try:
            # Set max to 100 aa (300 bp CDS)
            stats = filter_and_write_gff3(
                temp_gff3_file,
                output_file,
                min_protein_length=30,
                max_protein_length=100,
            )

            # Should keep genes 1, 2 (2 genes)
            # Should filter genes 3 (too small), 4 (too large), 5 (too large)
            assert stats["kept"] == 2
            assert stats["filtered"] == 3

            # Verify output file
            with open(output_file, "r") as f:
                content = f.read()
                assert "gene1" in content
                assert "gene2" in content
                assert "gene3" not in content
                assert "gene4" not in content
                assert "gene5" not in content
        finally:
            if os.path.exists(output_file):
                os.unlink(output_file)

    def test_filter_preserves_features(self, temp_gff3_file, temp_genome_fasta):
        """Test that filtering preserves all associated features."""
        output_file = tempfile.mktemp(suffix=".gff3")

        try:
            stats = filter_and_write_gff3(temp_gff3_file, output_file)

            # Verify that for kept genes, all features are present
            with open(output_file, "r") as f:
                lines = f.readlines()

                # Count feature types for gene1
                gene1_lines = [
                    l
                    for l in lines
                    if "gene1" in l or "mRNA1" in l or "CDS1" in l or "exon1" in l
                ]
                assert len(gene1_lines) == 4  # gene, mRNA, CDS, exon

                # Verify feature types
                feature_types = [l.split("\t")[2] for l in gene1_lines]
                assert "gene" in feature_types
                assert "mRNA" in feature_types
                assert "CDS" in feature_types
                assert "exon" in feature_types
        finally:
            if os.path.exists(output_file):
                os.unlink(output_file)

    def test_filter_empty_input(self, temp_genome_fasta):
        """Test filtering with an empty GFF3 file."""
        input_file = tempfile.mktemp(suffix=".gff3")
        output_file = tempfile.mktemp(suffix=".gff3")

        # Create empty input file
        open(input_file, "w").close()

        try:
            stats = filter_and_write_gff3(input_file, output_file)

            assert stats["kept"] == 0
            assert stats["filtered"] == 0
        finally:
            if os.path.exists(input_file):
                os.unlink(input_file)
            if os.path.exists(output_file):
                os.unlink(output_file)
