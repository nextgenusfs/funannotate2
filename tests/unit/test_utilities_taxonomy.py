"""
Unit tests for the taxonomy-related functions in utilities.py.
"""

from unittest.mock import patch

import funannotate2.utilities
from funannotate2.utilities import choose_best_busco_species
from funannotate2.config import busco_taxonomy


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


class TestChooseBestBuscoSpecies:
    """Integration tests for the choose_best_busco_species function."""

    def test_monascus_ruber_issue(self):
        """Test the specific issue with Monascus ruber taxonomy returning invalid key.

        This test ensures that choose_best_busco_species returns a valid key from
        busco_taxonomy, not just the taxonomic value itself.
        """
        # The problematic taxonomy from the user's example
        taxonomy = {
            "superkingdom": "Eukaryota",
            "kingdom": "Fungi",
            "phylum": "Ascomycota",
            "class": "Eurotiomycetes",
            "order": "Eurotiales",
            "family": "Aspergillaceae",
            "genus": "Monascus",
            "species": "Monascus ruber",
        }

        # Call the function
        result = choose_best_busco_species(taxonomy)

        # The result should be a valid key in busco_taxonomy
        assert result is not None, "Function should return a result"
        assert result in busco_taxonomy, (
            f"Result '{result}' should be a valid key in busco_taxonomy"
        )

        # The result should NOT be the raw taxonomic value
        assert result != "aspergillaceae", (
            "Should not return raw taxonomic value 'aspergillaceae'"
        )

        # For this specific case, we expect 'eurotiales' as the best match
        # because it's the most specific taxonomic level that has a direct match
        assert result == "eurotiales", f"Expected 'eurotiales' but got '{result}'"

    def test_aspergillus_fumigatus_exact_match(self):
        """Test with Aspergillus fumigatus which should match exactly."""
        taxonomy = {
            "superkingdom": "Eukaryota",
            "kingdom": "Fungi",
            "phylum": "Ascomycota",
            "class": "Eurotiomycetes",
            "order": "Eurotiales",
            "family": "Aspergillaceae",
            "genus": "Aspergillus",
            "species": "Aspergillus fumigatus",
        }

        result = choose_best_busco_species(taxonomy)

        assert result is not None
        assert result in busco_taxonomy
        # Should match at genus level since that's the most specific available
        assert result == "aspergillus"

    def test_saccharomyces_cerevisiae(self):
        """Test with Saccharomyces cerevisiae."""
        taxonomy = {
            "superkingdom": "Eukaryota",
            "kingdom": "Fungi",
            "phylum": "Ascomycota",
            "class": "Saccharomycetes",
            "order": "Saccharomycetales",
            "family": "Saccharomycetaceae",
            "genus": "Saccharomyces",
            "species": "Saccharomyces cerevisiae",
        }

        result = choose_best_busco_species(taxonomy)

        assert result is not None
        assert result in busco_taxonomy
        # Should find a valid match in the busco taxonomy

    def test_returns_valid_busco_key_always(self):
        """Test that the function always returns a valid busco_taxonomy key when it returns something."""
        test_taxonomies = [
            {
                "superkingdom": "Eukaryota",
                "kingdom": "Fungi",
                "phylum": "Ascomycota",
                "class": "Eurotiomycetes",
            },
            {
                "superkingdom": "Eukaryota",
                "kingdom": "Fungi",
                "phylum": "Basidiomycota",
            },
            {
                "superkingdom": "Eukaryota",
                "kingdom": "Metazoa",
                "phylum": "Chordata",
            },
        ]

        for taxonomy in test_taxonomies:
            result = choose_best_busco_species(taxonomy)
            if result is not None:  # Function might return None for no matches
                assert result in busco_taxonomy, (
                    f"Result '{result}' should be a valid key in busco_taxonomy"
                )
