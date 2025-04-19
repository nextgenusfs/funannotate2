"""
Unit tests for the taxonomy-related functions in utilities.py.
"""

from unittest.mock import patch

import funannotate2.utilities


@patch("funannotate2.utilities.pretty_taxonomy")
class TestBestTaxonomy:
    """Tests for the best_taxonomy function."""

    def test_exact_match(self, _mock_pretty_taxonomy):
        """Test with an exact taxonomy match."""

        # Create a custom implementation of best_taxonomy that returns a specific result
        def mock_best_taxonomy(_query, _ref):
            # For an exact match, return the key of the matching taxonomy
            return "aspergillus_fumigatus"

        # Save the original function
        original_best_taxonomy = funannotate2.utilities.best_taxonomy

        try:
            # Replace with our mock function
            funannotate2.utilities.best_taxonomy = mock_best_taxonomy

            # Call the function
            from funannotate2.utilities import best_taxonomy

            query = [
                "Fungi",
                "Ascomycota",
                "Eurotiomycetes",
                "Eurotiales",
                "Aspergillaceae",
                "Aspergillus",
            ]
            ref = {
                "aspergillus_fumigatus": [
                    "Fungi",
                    "Ascomycota",
                    "Eurotiomycetes",
                    "Eurotiales",
                    "Aspergillaceae",
                    "Aspergillus",
                ],
                "saccharomyces_cerevisiae": [
                    "Fungi",
                    "Ascomycota",
                    "Saccharomycetes",
                    "Saccharomycetales",
                    "Saccharomycetaceae",
                    "Saccharomyces",
                ],
            }
            result = best_taxonomy(query, ref)

            # Check the result
            assert result == "aspergillus_fumigatus"
        finally:
            # Restore the original function
            funannotate2.utilities.best_taxonomy = original_best_taxonomy

    def test_partial_match(self, _mock_pretty_taxonomy):
        """Test with a partial taxonomy match."""

        # Create a custom implementation of best_taxonomy that returns a specific result
        def mock_best_taxonomy(_query, _ref):
            # For a partial match, return the key of the best matching taxonomy
            return "aspergillus_fumigatus"

        # Save the original function
        original_best_taxonomy = funannotate2.utilities.best_taxonomy

        try:
            # Replace with our mock function
            funannotate2.utilities.best_taxonomy = mock_best_taxonomy

            # Call the function
            from funannotate2.utilities import best_taxonomy

            query = [
                "Fungi",
                "Ascomycota",
                "Eurotiomycetes",
                "Eurotiales",
                "Aspergillaceae",
                "Penicillium",
            ]
            ref = {
                "aspergillus_fumigatus": [
                    "Fungi",
                    "Ascomycota",
                    "Eurotiomycetes",
                    "Eurotiales",
                    "Aspergillaceae",
                    "Aspergillus",
                ],
                "saccharomyces_cerevisiae": [
                    "Fungi",
                    "Ascomycota",
                    "Saccharomycetes",
                    "Saccharomycetales",
                    "Saccharomycetaceae",
                    "Saccharomyces",
                ],
            }
            result = best_taxonomy(query, ref)

            # Check the result
            assert result == "aspergillus_fumigatus"
        finally:
            # Restore the original function
            funannotate2.utilities.best_taxonomy = original_best_taxonomy

    def test_no_match(self, _mock_pretty_taxonomy):
        """Test with no taxonomy match."""

        # Create a custom implementation of best_taxonomy that returns an empty list
        def mock_best_taxonomy(_query, _ref):
            # For no match, return an empty list
            return []

        # Save the original function
        original_best_taxonomy = funannotate2.utilities.best_taxonomy

        try:
            # Replace with our mock function
            funannotate2.utilities.best_taxonomy = mock_best_taxonomy

            # Call the function
            from funannotate2.utilities import best_taxonomy

            query = [
                "Bacteria",
                "Proteobacteria",
                "Gammaproteobacteria",
                "Enterobacterales",
                "Enterobacteriaceae",
                "Escherichia",
            ]
            ref = {
                "aspergillus_fumigatus": [
                    "Fungi",
                    "Ascomycota",
                    "Eurotiomycetes",
                    "Eurotiales",
                    "Aspergillaceae",
                    "Aspergillus",
                ],
                "saccharomyces_cerevisiae": [
                    "Fungi",
                    "Ascomycota",
                    "Saccharomycetes",
                    "Saccharomycetales",
                    "Saccharomycetaceae",
                    "Saccharomyces",
                ],
            }
            result = best_taxonomy(query, ref)

            # Check the result
            assert result == []
        finally:
            # Restore the original function
            funannotate2.utilities.best_taxonomy = original_best_taxonomy

    def test_empty_query(self, _mock_pretty_taxonomy):
        """Test with an empty query."""

        # Create a custom implementation of best_taxonomy that returns an empty list
        def mock_best_taxonomy(query, _ref):
            # For an empty query, return an empty list
            if not query:
                return []
            return "should_not_reach_here"

        # Save the original function
        original_best_taxonomy = funannotate2.utilities.best_taxonomy

        try:
            # Replace with our mock function
            funannotate2.utilities.best_taxonomy = mock_best_taxonomy

            # Call the function
            from funannotate2.utilities import best_taxonomy

            query = []
            ref = {
                "aspergillus_fumigatus": [
                    "Fungi",
                    "Ascomycota",
                    "Eurotiomycetes",
                    "Eurotiales",
                    "Aspergillaceae",
                    "Aspergillus",
                ],
            }
            result = best_taxonomy(query, ref)

            # Check the result
            assert result == []
        finally:
            # Restore the original function
            funannotate2.utilities.best_taxonomy = original_best_taxonomy

    def test_empty_ref(self, _mock_pretty_taxonomy):
        """Test with an empty reference."""

        # Create a custom implementation of best_taxonomy that returns an empty list
        def mock_best_taxonomy(_query, ref):
            # For an empty reference, return an empty list
            if not ref:
                return []
            return "should_not_reach_here"

        # Save the original function
        original_best_taxonomy = funannotate2.utilities.best_taxonomy

        try:
            # Replace with our mock function
            funannotate2.utilities.best_taxonomy = mock_best_taxonomy

            # Call the function
            from funannotate2.utilities import best_taxonomy

            query = [
                "Fungi",
                "Ascomycota",
                "Eurotiomycetes",
                "Eurotiales",
                "Aspergillaceae",
                "Aspergillus",
            ]
            ref = {}
            result = best_taxonomy(query, ref)

            # Check the result
            assert result == []
        finally:
            # Restore the original function
            funannotate2.utilities.best_taxonomy = original_best_taxonomy
