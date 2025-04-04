"""
Additional unit tests for the fastx module.
"""

import os
import pytest
from funannotate2.fastx import (
    simplify_headers,
    simplify_headers_drop,
    list2groups,
    fasta2chunks,
)


class TestSimplifyHeaders:
    """Tests for the simplify_headers function."""

    def test_simplify_headers(self, temp_dir, test_data_dir):
        """Test simplifying headers in a FASTA file."""
        # Setup input and output files
        fasta_file = os.path.join(test_data_dir, "test_headers.fasta")
        output_file = os.path.join(temp_dir, "simplified.fasta")

        # Run the function
        names = simplify_headers(fasta_file, output_file, base="contig_")

        # Check that the output file exists
        assert os.path.exists(output_file)

        # Check the mapping dictionary
        assert len(names) == 3
        assert "contig_1" in names
        assert "contig_2" in names
        assert "contig_3" in names
        # pyfastx only returns the part of the header before the first space
        assert names["contig_1"] == "scaffold_1"
        assert names["contig_2"] == "scaffold_2"
        assert names["contig_3"] == "scaffold_3"

        # Check the content of the output file
        with open(output_file, "r") as f:
            content = f.read()

        # Verify the headers are simplified
        assert ">contig_1" in content
        assert ">contig_2" in content
        assert ">contig_3" in content

        # Verify the original headers are not present
        assert ">scaffold_1" not in content
        assert ">scaffold_2" not in content
        assert ">scaffold_3" not in content


class TestSimplifyHeadersDrop:
    """Tests for the simplify_headers_drop function."""

    def test_simplify_headers_drop(self, temp_dir, test_data_dir):
        """Test simplifying headers and dropping specified contigs."""
        # Setup input and output files
        fasta_file = os.path.join(test_data_dir, "test_headers.fasta")
        keep_file = os.path.join(temp_dir, "keep.fasta")
        drop_file = os.path.join(temp_dir, "drop.fasta")

        # Run the function with scaffold_2 to be dropped
        names = simplify_headers_drop(
            fasta_file,
            keep_file,
            drop_file,
            base="contig_",
            drop=["scaffold_2"],
        )

        # Check that the output files exist
        assert os.path.exists(keep_file)
        assert os.path.exists(drop_file)

        # Check the mapping dictionary - note that the drop list only contains "scaffold_2",
        # but pyfastx only returns the part of the header before the first space
        assert len(names) == 2
        assert "contig_1" in names
        assert "contig_3" in names
        assert names["contig_1"] == "scaffold_1"
        assert names["contig_3"] == "scaffold_3"

        # Check the content of the keep file
        with open(keep_file, "r") as f:
            keep_content = f.read()

        # Verify the headers are simplified and scaffold_2 is not present
        assert ">contig_1" in keep_content
        assert ">contig_3" in keep_content
        assert ">scaffold_2" not in keep_content
        assert "contig_2" not in keep_content

        # Check the content of the drop file
        with open(drop_file, "r") as f:
            drop_content = f.read()

        # Verify scaffold_2 is in the drop file with its original header
        assert ">scaffold_2" in drop_content


class TestList2Groups:
    """Tests for the list2groups function."""

    def test_empty_list(self):
        """Test with an empty list."""
        result = list(list2groups([]))
        assert result == []

    def test_single_number(self):
        """Test with a single number."""
        result = list(list2groups([5]))
        assert result == [(5, 5)]

    def test_continuous_numbers(self):
        """Test with continuous numbers."""
        result = list(list2groups([1, 2, 3, 4, 5]))
        assert result == [(1, 5)]

    def test_multiple_groups(self):
        """Test with multiple groups of continuous numbers."""
        result = list(list2groups([1, 2, 3, 5, 6, 8, 9, 10]))
        assert result == [(1, 3), (5, 6), (8, 10)]

    def test_unsorted_numbers(self):
        """Test with unsorted numbers."""
        result = list(list2groups([3, 1, 2, 10, 8, 9, 5, 6]))
        # list2groups assumes the input is sorted
        assert result != [(1, 3), (5, 6), (8, 10)]


class TestFasta2Chunks:
    """Tests for the fasta2chunks function."""

    def test_fasta2chunks(self, temp_dir, test_data_dir):
        """Test splitting a FASTA file into chunks."""
        # Setup input file
        fasta_file = os.path.join(test_data_dir, "test_headers.fasta")
        output_dir = os.path.join(temp_dir, "chunks")

        # Run the function to split into 2 chunks
        chunk_files = fasta2chunks(fasta_file, 2, output_dir)

        # Check that the output directory exists
        assert os.path.isdir(output_dir)

        # Check that the correct number of chunk files were created
        assert len(chunk_files) == 2

        # Check that all chunk files exist
        for chunk_file in chunk_files:
            assert os.path.exists(chunk_file)

        # Check the content of the first chunk file
        with open(chunk_files[0], "r") as f:
            chunk1_content = f.read()

        # The first chunk should contain at least one sequence
        assert ">" in chunk1_content

        # Check that the second chunk file exists
        # The second chunk might be empty if there are only 2 sequences,
        # but the file should exist
        assert os.path.exists(chunk_files[1])
