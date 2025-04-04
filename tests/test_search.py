#!/usr/bin/env python3

import os
import unittest
import tempfile
from unittest.mock import patch, MagicMock
from buscolite.search import (
    tblastn_version,
    miniprot_version,
    pyhmmer_version,
    merge_overlapping_hits,
)


class TestSearchFunctions(unittest.TestCase):
    """Test basic functions in the search module."""

    @patch("subprocess.Popen")
    def test_tblastn_version(self, mock_popen):
        """Test getting tblastn version."""
        # Mock the subprocess.Popen to return a specific version
        process_mock = MagicMock()
        process_mock.communicate.return_value = ("tblastn: 2.10.1+", "")
        mock_popen.return_value = process_mock

        version = tblastn_version()
        self.assertEqual(version, "2.10.1")

        # Check that Popen was called with the correct arguments
        mock_popen.assert_called_once_with(
            ["tblastn", "-version"],
            stdout=-1,
            stderr=-1,
            universal_newlines=True,
        )

    @patch("subprocess.Popen")
    def test_miniprot_version(self, mock_popen):
        """Test getting miniprot version."""
        # Mock the subprocess.Popen to return a specific version
        process_mock = MagicMock()
        process_mock.communicate.return_value = ("0.7", "")
        mock_popen.return_value = process_mock

        version = miniprot_version()
        self.assertEqual(version, "0.7")

        # Check that Popen was called with the correct arguments
        mock_popen.assert_called_once_with(
            ["miniprot", "--version"],
            stdout=-1,
            stderr=-1,
            universal_newlines=True,
        )

    @patch("buscolite.search.pyhmmer.__version__", "0.10.15")
    def test_pyhmmer_version(self):
        """Test getting pyhmmer version."""
        version = pyhmmer_version()
        self.assertEqual(version, "0.10.15")

    def test_merge_overlapping_hits_single(self):
        """Test merging overlapping hits with a single hit."""
        hits = [{"coords": (100, 200), "score": 10}]
        result = merge_overlapping_hits(hits)
        self.assertEqual(result, hits)

    def test_merge_overlapping_hits_overlapping(self):
        """Test merging overlapping hits."""
        hits = [
            {"coords": (100, 200), "score": 10},
            {"coords": (150, 250), "score": 20},
            {"coords": (300, 400), "score": 30},
        ]
        result = merge_overlapping_hits(hits)
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0]["coords"], (100, 250))
        self.assertEqual(result[1]["coords"], (300, 400))

    def test_merge_overlapping_hits_nearby(self):
        """Test merging nearby hits within fluff distance."""
        hits = [
            {"coords": (100, 200), "score": 10},
            {"coords": (5000, 6000), "score": 20},
            {"coords": (210, 300), "score": 30},  # Within default fluff (10000)
        ]
        result = merge_overlapping_hits(hits)
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0]["coords"], (100, 300))
        self.assertEqual(result[1]["coords"], (5000, 6000))

    def test_merge_overlapping_hits_custom_fluff(self):
        """Test merging nearby hits with custom fluff distance."""
        hits = [
            {"coords": (100, 200), "score": 10},
            {"coords": (250, 350), "score": 20},  # 50 bp away from previous
        ]
        # With fluff=100, they should merge
        result = merge_overlapping_hits(hits, fluff=100)
        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]["coords"], (100, 350))

        # With fluff=20, they should not merge
        result = merge_overlapping_hits(hits, fluff=20)
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0]["coords"], (100, 200))
        self.assertEqual(result[1]["coords"], (250, 350))


if __name__ == "__main__":
    unittest.main()
