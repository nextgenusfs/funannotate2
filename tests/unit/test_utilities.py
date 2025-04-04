"""
Unit tests for the utilities module.
"""
import os
import pytest
from funannotate2.utilities import (
    merge_coordinates,
    naming_slug,
    create_tmpdir,
    readBlocks,
)


class TestMergeCoordinates:
    """Tests for the merge_coordinates function."""

    def test_empty_list(self):
        """Test with an empty list."""
        assert merge_coordinates([]) == []

    def test_single_interval(self):
        """Test with a single interval."""
        assert merge_coordinates([[1, 3]]) == [[1, 3]]

    def test_non_overlapping_intervals(self):
        """Test with non-overlapping intervals."""
        assert merge_coordinates([[1, 3], [5, 7], [9, 11]]) == [[1, 3], [5, 7], [9, 11]]

    def test_overlapping_intervals(self):
        """Test with overlapping intervals."""
        assert merge_coordinates([[1, 3], [2, 6], [8, 10], [15, 18]]) == [
            [1, 6],
            [8, 10],
            [15, 18],
        ]

    def test_completely_overlapping_intervals(self):
        """Test with completely overlapping intervals."""
        assert merge_coordinates([[1, 10], [2, 5], [3, 7]]) == [[1, 10]]

    def test_adjacent_intervals(self):
        """Test with adjacent intervals."""
        assert merge_coordinates([[1, 3], [3, 5], [7, 9]]) == [[1, 5], [7, 9]]

    def test_unsorted_intervals(self):
        """Test with unsorted intervals."""
        assert merge_coordinates([[8, 10], [1, 3], [2, 6], [15, 18]]) == [
            [1, 6],
            [8, 10],
            [15, 18],
        ]


class TestNamingSlug:
    """Tests for the naming_slug function."""

    def test_basic_slug(self):
        """Test basic slug generation."""
        assert (
            naming_slug("aspergillus fumigatus", "Af293")
            == "Aspergillus_fumigatus_Af293"
        )

    def test_lowercase_slug(self):
        """Test lowercase slug generation."""
        assert (
            naming_slug("Aspergillus fumigatus", "Af293", lowercase=True)
            == "aspergillus_fumigatus_Af293"
        )

    def test_no_strain(self):
        """Test slug generation without strain."""
        assert naming_slug("Aspergillus fumigatus", None) == "Aspergillus_fumigatus"

    def test_strain_with_spaces(self):
        """Test slug generation with spaces in strain."""
        assert (
            naming_slug("Aspergillus fumigatus", "Af 293")
            == "Aspergillus_fumigatus_Af293"
        )


class TestCreateTmpdir:
    """Tests for the create_tmpdir function."""

    def test_create_tmpdir_with_outdir(self, temp_dir):
        """Test creating a temporary directory with an output directory."""
        tmpdir = create_tmpdir(temp_dir, base="test")
        assert os.path.exists(tmpdir)
        assert os.path.isdir(tmpdir)
        assert tmpdir.startswith(os.path.abspath(temp_dir))

    def test_create_tmpdir_with_tmp(self):
        """Test creating a temporary directory in /tmp."""
        tmpdir = create_tmpdir("/tmp", base="test")
        assert os.path.exists(tmpdir)
        assert os.path.isdir(tmpdir)
        assert tmpdir.startswith("/tmp")

    def test_create_tmpdir_without_outdir(self):
        """Test creating a temporary directory without an output directory."""
        tmpdir = create_tmpdir(None, base="test")
        assert os.path.exists(tmpdir)
        assert os.path.isdir(tmpdir)
        # Clean up
        os.rmdir(tmpdir)


class TestReadBlocks:
    """Tests for the readBlocks function."""

    def test_read_blocks(self):
        """Test reading blocks from a source."""
        source = [
            "# Block 1",
            "Line 1",
            "Line 2",
            "# Block 2",
            "Line 3",
            "Line 4",
            "# Block 3",
            "Line 5",
        ]
        blocks = list(readBlocks(source, "#"))
        assert len(blocks) == 3
        assert blocks[0] == ["# Block 1", "Line 1", "Line 2"]
        assert blocks[1] == ["# Block 2", "Line 3", "Line 4"]
        assert blocks[2] == ["# Block 3", "Line 5"]

    def test_read_blocks_empty_source(self):
        """Test reading blocks from an empty source."""
        blocks = list(readBlocks([], "#"))
        assert len(blocks) == 1
        assert blocks[0] == []

    def test_read_blocks_no_pattern(self):
        """Test reading blocks with no pattern match."""
        source = ["Line 1", "Line 2", "Line 3"]
        blocks = list(readBlocks(source, "#"))
        assert len(blocks) == 1
        assert blocks[0] == ["Line 1", "Line 2", "Line 3"]
