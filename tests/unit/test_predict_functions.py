"""
Unit tests for specific functions in the predict module.
"""

from unittest.mock import patch

import pytest

import funannotate2.predict


class TestSanitizeExternal:
    """Tests for the sanitize_external function."""

    @pytest.fixture
    def mock_gff2dict(self):
        """Mock the gff2dict function."""
        with patch("funannotate2.predict.gff2dict") as mock:
            # Set up the mock to return a sample gene dictionary
            mock.return_value = {
                "gene1": {
                    "contig": "original_contig1",
                    "source": "Augustus",
                    "location": [1, 1000],
                    "strand": "+",
                    "type": "gene",
                },
                "gene2": {
                    "contig": "original_contig2",
                    "source": "GeneMark",
                    "location": [2000, 3000],
                    "strand": "-",
                    "type": "gene",
                },
            }
            yield mock

    def test_sanitize_external(self, mock_gff2dict):
        """Test sanitizing external gene annotations."""
        # Set up test data
        gff3_file = "test.gff3"
        fasta_file = "test.fasta"
        contig_map = {
            "new_contig1": "original_contig1",
            "new_contig2": "original_contig2",
        }

        # Call the function
        result, sources = funannotate2.predict.sanitize_external(gff3_file, fasta_file, contig_map)

        # Check that gff2dict was called with the correct arguments
        mock_gff2dict.assert_called_once_with(gff3_file, fasta_file, debug=False)

        # Check the result
        assert len(result) == 2
        assert result["gene1"]["contig"] == "new_contig1"
        assert result["gene2"]["contig"] == "new_contig2"
        assert result["gene1"]["source"] == "augustus"
        assert result["gene2"]["source"] == "genemark"

        # Check the sources
        assert sources == {"augustus", "genemark"}

    def test_sanitize_external_with_missing_contig(self, mock_gff2dict):
        """Test sanitizing external gene annotations with a missing contig in the map."""
        # Set up test data
        gff3_file = "test.gff3"
        fasta_file = "test.fasta"
        contig_map = {
            "new_contig1": "original_contig1",
            # Missing mapping for original_contig2
        }

        # Call the function
        result, sources = funannotate2.predict.sanitize_external(gff3_file, fasta_file, contig_map)

        # Check the result
        assert len(result) == 2
        assert result["gene1"]["contig"] == "new_contig1"
        assert result["gene2"]["contig"] is None  # No mapping for original_contig2


class TestMergeRenameModels:
    """Tests for the merge_rename_models function."""

    @pytest.fixture
    def mock_gff2dict(self):
        """Mock the gff2dict function."""
        with patch("funannotate2.predict.gff2dict") as mock:
            # Set up the mock to return sample gene dictionaries for different GFF files
            mock.side_effect = lambda gff, genome, annotation=None: {
                **(annotation or {}),
                f"gene_{gff}": {
                    "contig": f"contig_{gff}",
                    "location": [1, 1000],
                    "strand": "+",
                    "type": "gene",
                    "ids": [f"gene_{gff}-T1"],
                },
            }
            yield mock

    @pytest.fixture
    def mock_dict2gff3(self):
        """Mock the dict2gff3 function."""
        with patch("funannotate2.predict.dict2gff3") as mock:
            yield mock

    @pytest.fixture
    def mock_natsorted(self):
        """Mock the natsorted function."""
        with patch("funannotate2.predict.natsorted") as mock:
            # Just return the input for simplicity
            mock.side_effect = lambda x, key: sorted(x, key=key)
            yield mock

    def test_merge_rename_models(self, mock_gff2dict, mock_dict2gff3, mock_natsorted):
        """Test merging and renaming gene models."""
        # Set up test data
        gff_list = ["gff1", "gff2"]
        genome = "genome.fasta"
        output = "output.gff3"
        locus_tag = "TEST_"
        contig_map = {"new_contig1": "contig_gff1", "new_contig2": "contig_gff2"}

        # Call the function
        result = funannotate2.predict.merge_rename_models(
            gff_list, genome, output, locus_tag, contig_map
        )

        # Check that gff2dict was called for each GFF file
        assert mock_gff2dict.call_count == 2

        # Check that dict2gff3 was called with the correct arguments
        mock_dict2gff3.assert_called_once()
        args, kwargs = mock_dict2gff3.call_args
        assert kwargs["output"] == output
        assert kwargs["source"] == "funannotate2"

        # Check the result
        assert len(result) == 2
        # Since we're using a mock implementation, we can't check the exact values
        # Just check that the keys are in the expected format
        assert any(key.startswith("TEST_") for key in result.keys())

    def test_merge_rename_models_with_empty_list(
        self, mock_gff2dict, mock_dict2gff3, mock_natsorted
    ):
        """Test merging and renaming gene models with an empty GFF list."""
        # Set up test data
        gff_list = []
        genome = "genome.fasta"
        output = "output.gff3"

        # Call the function
        result = funannotate2.predict.merge_rename_models(gff_list, genome, output)

        # Check that gff2dict was not called
        mock_gff2dict.assert_not_called()

        # Check that dict2gff3 was called with an empty dictionary
        mock_dict2gff3.assert_called_once()
        args, kwargs = mock_dict2gff3.call_args
        assert kwargs["output"] == output
        assert kwargs["source"] == "funannotate2"

        # Check the result
        assert len(result) == 0
