"""
Unit tests for the fastx module.
"""

import os

from funannotate2.fastx import (
    annotate_fasta,
    countfasta,
    fasta2dict,
    mergefasta,
    softwrap,
)


class TestSoftwrap:
    """Tests for the softwrap function."""

    def test_empty_string(self):
        """Test with an empty string."""
        assert softwrap("") == ""

    def test_short_string(self):
        """Test with a string shorter than the wrap length."""
        assert softwrap("ATGC") == "ATGC"

    def test_exact_length_string(self):
        """Test with a string exactly the wrap length."""
        assert softwrap("A" * 80) == "A" * 80

    def test_long_string(self):
        """Test with a string longer than the wrap length."""
        assert softwrap("A" * 100) == "A" * 80 + "\n" + "A" * 20

    def test_custom_wrap_length(self):
        """Test with a custom wrap length."""
        assert softwrap("A" * 100, every=50) == "A" * 50 + "\n" + "A" * 50


class TestCountfasta:
    """Tests for the countfasta function."""

    def test_count_fasta(self, test_data_dir):
        """Test counting sequences in a FASTA file."""
        fasta_file = os.path.join(test_data_dir, "test.fasta")
        assert countfasta(fasta_file) == 3


class TestFasta2dict:
    """Tests for the fasta2dict function."""

    def test_fasta2dict(self, test_data_dir):
        """Test converting a FASTA file to a dictionary."""
        fasta_file = os.path.join(test_data_dir, "test.fasta")
        fasta_dict = fasta2dict(fasta_file)

        # Check that we have the expected number of sequences
        assert len(fasta_dict) == 3

        # Check that the keys are correct
        assert "contig1" in fasta_dict
        assert "contig2" in fasta_dict
        assert "contig3" in fasta_dict

        # Check that the sequences are correct (just check the first 10 characters)
        assert fasta_dict["contig1"][:10] == "ATGCGTACGT"
        assert fasta_dict["contig2"][:10] == "ATGCGTACGT"
        assert fasta_dict["contig3"][:10] == "ATGCGTACGT"


class TestAnnotateFasta:
    """Tests for the annotate_fasta function."""

    def test_annotate_fasta(self, temp_dir, test_data_dir):
        """Test annotating sequences in a FASTA file."""
        # Setup input and output files
        fasta_file = os.path.join(test_data_dir, "test.fasta")
        output_file = os.path.join(temp_dir, "annotated.fasta")

        # Run the function with specific IDs to annotate
        annotate_fasta(fasta_file, output_file, ids=["contig1", "contig3"])

        # Check that the output file exists
        assert os.path.exists(output_file)

        # Read the output file directly to check the headers
        with open(output_file, "r") as f:
            content = f.read()

        # Check that the annotations are present in the file content
        assert ">contig1 [mcode=4] [location=mitochondrion]" in content
        assert ">contig3 [mcode=4] [location=mitochondrion]" in content
        assert ">contig2\n" in content  # This one should not be annotated

    def test_annotate_fasta_custom_annotation(self, temp_dir, test_data_dir):
        """Test annotating sequences with a custom annotation."""
        # Setup input and output files
        fasta_file = os.path.join(test_data_dir, "test.fasta")
        output_file = os.path.join(temp_dir, "custom_annotated.fasta")

        # Run the function with a custom annotation
        custom_annotation = "[custom=test]"
        annotate_fasta(fasta_file, output_file, ids=["contig2"], annotation=custom_annotation)

        # Check that the output file exists
        assert os.path.exists(output_file)

        # Read the output file directly to check the headers
        with open(output_file, "r") as f:
            content = f.read()

        # Check that the annotations are present in the file content
        assert ">contig2 [custom=test]" in content
        assert ">contig1\n" in content  # This one should not be annotated
        assert ">contig3\n" in content  # This one should not be annotated


class TestMergefasta:
    """Tests for the mergefasta function."""

    def test_mergefasta(self, temp_dir, test_data_dir):
        """Test merging and dereplicating FASTA files."""
        # Setup input and output files
        fasta_file1 = os.path.join(test_data_dir, "test1.fasta")
        fasta_file2 = os.path.join(test_data_dir, "test2.fasta")
        output_file = os.path.join(temp_dir, "merged.fasta")

        # Run the function to merge the files
        _, total_count = mergefasta([fasta_file1, fasta_file2], output_file)

        # Check that the output file exists
        assert os.path.exists(output_file)

        # Check the counts - there are 5 total sequences in the input files
        assert total_count == 5

        # Read the output file to check its contents
        with open(output_file, "r") as f:
            content = f.read()

        # Check that the expected sequences are in the output
        # The function deduplicates based on sequence content, not ID
        # seq1, seq2, and seq3 have identical sequences, so only one is kept
        assert ">seq1" in content  # seq1 is kept as it's the first occurrence
        assert ">seq4" in content  # seq4 has a unique sequence

        # seq2 and seq3 should be deduplicated and not present
        assert ">seq2" not in content
        assert ">seq3" not in content
