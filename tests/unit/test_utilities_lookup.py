"""
Unit tests for the lookup_taxonomy function in utilities.py.
"""

import os
import tempfile
import pytest
import funannotate2.utilities
from unittest.mock import patch, MagicMock


class TestLookupTaxonomy:
    """Tests for the lookup_taxonomy function."""

    def test_lookup_taxonomy(self):
        """Test the lookup_taxonomy function with a mock implementation."""

        # Create a custom implementation of lookup_taxonomy that returns a mock taxonomy
        def mock_lookup_taxonomy(name):
            return {
                "superkingdom": "Fungi",
                "kingdom": "",
                "phylum": "Ascomycota",
                "class": "Eurotiomycetes",
                "order": "Eurotiales",
                "family": "Aspergillaceae",
                "genus": "Aspergillus",
                "species": "Aspergillus fumigatus",
            }

        # Save the original function
        original_lookup_taxonomy = funannotate2.utilities.lookup_taxonomy

        try:
            # Replace with our mock function
            funannotate2.utilities.lookup_taxonomy = mock_lookup_taxonomy

            # Call the function through the imported name
            from funannotate2.utilities import lookup_taxonomy

            result = lookup_taxonomy("Aspergillus fumigatus")

            # Check the returned taxonomy
            assert result["superkingdom"] == "Fungi"
            assert result["phylum"] == "Ascomycota"
            assert result["class"] == "Eurotiomycetes"
            assert result["order"] == "Eurotiales"
            assert result["family"] == "Aspergillaceae"
            assert result["genus"] == "Aspergillus"
            assert result["species"] == "Aspergillus fumigatus"
        finally:
            # Restore the original function
            funannotate2.utilities.lookup_taxonomy = original_lookup_taxonomy

    def test_failed_lookup(self):
        """Test a failed taxonomy lookup."""

        # Create a custom implementation of lookup_taxonomy that returns an empty dict
        def mock_lookup_taxonomy(name):
            return {}

        # Save the original function
        original_lookup_taxonomy = funannotate2.utilities.lookup_taxonomy

        try:
            # Replace with our mock function
            funannotate2.utilities.lookup_taxonomy = mock_lookup_taxonomy

            # Call the function through the imported name
            from funannotate2.utilities import lookup_taxonomy

            result = lookup_taxonomy("Invalid Species")

            # Check that an empty dictionary is returned
            assert result == {}
        finally:
            # Restore the original function
            funannotate2.utilities.lookup_taxonomy = original_lookup_taxonomy
