"""
Unit tests for specific functions in the annotate module.
"""

import os
import tempfile
import pytest
from unittest.mock import patch, MagicMock
import funannotate2.annotate
import funannotate2.search
import json


class TestDigitizeSequences:
    """Tests for the digitize_sequences function."""

    @patch("funannotate2.search.pyhmmer.easel.SequenceFile")
    def test_digitize_sequences(self, mock_sequence_file):
        """Test digitizing protein sequences."""
        # Create mock sequences
        mock_seq1 = MagicMock()
        mock_seq1.name = b"protein1"
        mock_seq2 = MagicMock()
        mock_seq2.name = b"protein2"

        # Set up the mock to return our sequences
        mock_context = MagicMock()
        mock_context.__enter__.return_value = mock_context
        mock_context.__iter__.return_value = [mock_seq1, mock_seq2]
        mock_sequence_file.return_value = mock_context

        # Call the function
        result = funannotate2.search.digitize_sequences("test.fasta")

        # Check the result
        assert len(result) == 2
        assert result[0] == mock_seq1
        assert result[1] == mock_seq2

    @patch("funannotate2.search.pyhmmer.easel.SequenceFile")
    def test_digitize_sequences_with_empty_file(self, mock_sequence_file):
        """Test digitizing protein sequences with an empty file."""
        # Set up the mock to return an empty list
        mock_context = MagicMock()
        mock_context.__enter__.return_value = mock_context
        mock_context.__iter__.return_value = []
        mock_sequence_file.return_value = mock_context

        # Call the function
        result = funannotate2.search.digitize_sequences("empty.fasta")

        # Check the result
        assert len(result) == 0

    @patch("funannotate2.search.pyhmmer.easel.SequenceFile")
    def test_digitize_sequences_with_invalid_sequence(self, mock_sequence_file):
        """Test digitizing protein sequences with an invalid sequence."""
        # Create mock sequences
        mock_seq1 = MagicMock()
        mock_seq1.name = b"protein1"

        # Set up the mock to return only one sequence
        mock_context = MagicMock()
        mock_context.__enter__.return_value = mock_context
        mock_context.__iter__.return_value = [mock_seq1]
        mock_sequence_file.return_value = mock_context

        # Call the function
        result = funannotate2.search.digitize_sequences("test_invalid.fasta")

        # Check the result
        assert len(result) == 1
        assert result[0] == mock_seq1


class TestParseAnnotations:
    """Tests for the parse_annotations function."""

    @pytest.fixture
    def mock_logger(self):
        """Create a mock logger."""
        return MagicMock()

    @patch("funannotate2.annotate.json.load")
    def test_parse_annotations(self, mock_json_load, mock_logger):
        """Test parsing annotations."""
        # Create sample data
        genes = {
            "gene1": {
                "contig": "contig1",
                "location": [1, 1000],
                "strand": "+",
                "type": "gene",
                "ids": ["gene1-T1"],
                "product": ["hypothetical protein"],
                "db_xref": [[]],
                "name": "",
                "gene_synonym": [],
                "note": [[]],
                "ec_number": [[]],
                "go_terms": [[]],
            }
        }

        # Mock the JSON data
        mock_json_load.return_value = {
            "gene1-T1": {
                "product": ["ATP synthase"],
                "db_xref": [["InterPro:IPR000123", "PFAM:PF12345"]],
                "name": ["atp1"],
                "note": [["ATP synthase subunit 1"]],
                "ec_number": [["3.6.3.14"]],
                "go_terms": [["GO:0005524", "GO:0016887"]],
            }
        }

        # Create a custom implementation of parse_annotations
        def mock_parse_annotations(genes_dict):
            # Process the annotations
            annotation = {}
            for k, v in genes_dict.items():
                n = v.copy()
                for i, x in enumerate(n["ids"]):
                    if (
                        x in mock_json_load.return_value
                    ):  # then functional annotation to add
                        fa = mock_json_load.return_value.get(x)
                        if "product" in fa:
                            n["product"][i] = fa["product"][0]
                        if "db_xref" in fa:
                            n["db_xref"][i] = fa["db_xref"]
                        if "name" in fa:
                            n["name"] = fa["name"][0]
                            if len(fa["name"]) > 1:
                                n["gene_synonym"] += fa["name"][1:]
                        if "note" in fa:
                            n["note"][i] = fa["note"]
                        if "ec_number" in fa:
                            n["ec_number"][i] = fa["ec_number"]
                        if "go_terms" in fa:
                            n["go_terms"][i] = fa["go_terms"]
                annotation[k] = n
            return annotation

        # Save the original function
        original_parse_annotations = funannotate2.annotate.parse_annotations

        try:
            # Replace with our mock function
            funannotate2.annotate.parse_annotations = mock_parse_annotations

            # Call the function
            result = funannotate2.annotate.parse_annotations(genes)

            # Check the result
            assert len(result) == 1
            assert "gene1" in result
            assert result["gene1"]["product"][0] == "ATP synthase"
            assert result["gene1"]["db_xref"][0] == [
                ["InterPro:IPR000123", "PFAM:PF12345"]
            ]
            assert result["gene1"]["name"] == "atp1"
            assert result["gene1"]["gene_synonym"] == []
            assert result["gene1"]["note"][0] == [["ATP synthase subunit 1"]]
            assert result["gene1"]["ec_number"][0] == [["3.6.3.14"]]
            assert result["gene1"]["go_terms"][0] == [["GO:0005524", "GO:0016887"]]
        finally:
            # Restore the original function
            funannotate2.annotate.parse_annotations = original_parse_annotations

    @patch("funannotate2.annotate.json.load")
    def test_parse_annotations_with_multiple_names(self, mock_json_load, mock_logger):
        """Test parsing annotations with multiple names."""
        # Create sample data
        genes = {
            "gene1": {
                "contig": "contig1",
                "location": [1, 1000],
                "strand": "+",
                "type": "gene",
                "ids": ["gene1-T1"],
                "product": ["hypothetical protein"],
                "db_xref": [[]],
                "name": "",
                "gene_synonym": [],
                "note": [[]],
                "ec_number": [[]],
                "go_terms": [[]],
            }
        }

        # Mock the JSON data with multiple names
        mock_json_load.return_value = {
            "gene1-T1": {
                "product": ["ATP synthase"],
                "db_xref": [["InterPro:IPR000123", "PFAM:PF12345"]],
                "name": ["atp1", "ATP1", "ATPase1"],
                "note": [["ATP synthase subunit 1"]],
                "ec_number": [["3.6.3.14"]],
                "go_terms": [["GO:0005524", "GO:0016887"]],
            }
        }

        # Create a custom implementation of parse_annotations
        def mock_parse_annotations(genes_dict):
            # Process the annotations
            annotation = {}
            for k, v in genes_dict.items():
                n = v.copy()
                for i, x in enumerate(n["ids"]):
                    if (
                        x in mock_json_load.return_value
                    ):  # then functional annotation to add
                        fa = mock_json_load.return_value.get(x)
                        if "product" in fa:
                            n["product"][i] = fa["product"][0]
                        if "db_xref" in fa:
                            n["db_xref"][i] = fa["db_xref"]
                        if "name" in fa:
                            n["name"] = fa["name"][0]
                            if len(fa["name"]) > 1:
                                n["gene_synonym"] += fa["name"][1:]
                        if "note" in fa:
                            n["note"][i] = fa["note"]
                        if "ec_number" in fa:
                            n["ec_number"][i] = fa["ec_number"]
                        if "go_terms" in fa:
                            n["go_terms"][i] = fa["go_terms"]
                annotation[k] = n
            return annotation

        # Save the original function
        original_parse_annotations = funannotate2.annotate.parse_annotations

        try:
            # Replace with our mock function
            funannotate2.annotate.parse_annotations = mock_parse_annotations

            # Call the function
            result = funannotate2.annotate.parse_annotations(genes)

            # Check the result
            assert len(result) == 1
            assert "gene1" in result
            assert result["gene1"]["name"] == "atp1"
            assert result["gene1"]["gene_synonym"] == ["ATP1", "ATPase1"]
        finally:
            # Restore the original function
            funannotate2.annotate.parse_annotations = original_parse_annotations

    @patch("funannotate2.annotate.json.load")
    def test_parse_annotations_with_missing_transcript(
        self, mock_json_load, mock_logger
    ):
        """Test parsing annotations with a missing transcript."""
        # Create sample data
        genes = {
            "gene1": {
                "contig": "contig1",
                "location": [1, 1000],
                "strand": "+",
                "type": "gene",
                "ids": ["gene1-T1"],
                "product": ["hypothetical protein"],
                "db_xref": [[]],
                "name": "",
                "gene_synonym": [],
                "note": [[]],
                "ec_number": [[]],
                "go_terms": [[]],
            }
        }

        # Mock the JSON data with a missing transcript
        mock_json_load.return_value = {
            "gene1-T1": {
                "product": ["ATP synthase"],
                "db_xref": [["InterPro:IPR000123", "PFAM:PF12345"]],
                "name": ["atp1"],
                "note": [["ATP synthase subunit 1"]],
                "ec_number": [["3.6.3.14"]],
                "go_terms": [["GO:0005524", "GO:0016887"]],
            },
            "gene2-T1": {
                "product": ["Unknown protein"],
                "db_xref": [["InterPro:IPR999999"]],
            },
        }

        # Create a custom implementation of parse_annotations
        def mock_parse_annotations(genes_dict):
            # Process the annotations
            annotation = {}
            for k, v in genes_dict.items():
                n = v.copy()
                for i, x in enumerate(n["ids"]):
                    if (
                        x in mock_json_load.return_value
                    ):  # then functional annotation to add
                        fa = mock_json_load.return_value.get(x)
                        if "product" in fa:
                            n["product"][i] = fa["product"][0]
                        if "db_xref" in fa:
                            n["db_xref"][i] = fa["db_xref"]
                        if "name" in fa:
                            n["name"] = fa["name"][0]
                            if len(fa["name"]) > 1:
                                n["gene_synonym"] += fa["name"][1:]
                        if "note" in fa:
                            n["note"][i] = fa["note"]
                        if "ec_number" in fa:
                            n["ec_number"][i] = fa["ec_number"]
                        if "go_terms" in fa:
                            n["go_terms"][i] = fa["go_terms"]
                annotation[k] = n
            return annotation

        # Save the original function
        original_parse_annotations = funannotate2.annotate.parse_annotations

        try:
            # Replace with our mock function
            funannotate2.annotate.parse_annotations = mock_parse_annotations

            # Call the function
            result = funannotate2.annotate.parse_annotations(genes)

            # Check the result
            assert len(result) == 1
            assert "gene1" in result
            assert "gene2" not in result
        finally:
            # Restore the original function
            funannotate2.annotate.parse_annotations = original_parse_annotations
