"""
Integration tests for the buscolite Python API.
"""
import os
import sys
import json
import tempfile
import shutil
import pytest
from pathlib import Path
from buscolite.busco import (
    runbusco,
    check_lineage,
    load_config,
    load_cutoffs,
)


class TestBUSCOliteAPI:
    """Test the buscolite Python API."""

    @pytest.fixture(scope="class")
    def test_data_dir(self):
        """Return the path to the test data directory."""
        return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

    @pytest.fixture(scope="function")
    def temp_dir(self):
        """Create a temporary directory for test outputs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield tmpdir

    def test_check_lineage(self, test_data_dir):
        """Test the check_lineage function."""
        lineage_path = os.path.join(test_data_dir, "mock_lineage")
        valid, message = check_lineage(lineage_path)
        assert valid is True
        assert message == ""

    def test_load_config(self, test_data_dir):
        """Test the load_config function."""
        lineage_path = os.path.join(test_data_dir, "mock_lineage")
        config = load_config(lineage_path)
        assert config["name"] == "mock_lineage"
        assert config["species"] == "anidulans"
        assert config["domain"] == "eukaryota"
        assert config["creation_date"] == "2023-01-01"
        assert config["number_of_BUSCOs"] == "2"
        assert config["number_of_species"] == "1"

    def test_load_cutoffs(self, test_data_dir):
        """Test the load_cutoffs function."""
        lineage_path = os.path.join(test_data_dir, "mock_lineage")
        cutoffs = load_cutoffs(lineage_path)
        assert "busco1" in cutoffs
        assert "busco2" in cutoffs
        assert cutoffs["busco1"]["score"] == 100.0
        assert cutoffs["busco2"]["score"] == 200.0
        assert cutoffs["busco1"]["length"] == 300
        assert cutoffs["busco2"]["length"] == 400
        assert cutoffs["busco1"]["sigma"] == 1.5
        assert cutoffs["busco2"]["sigma"] == 1.0  # Default value when 0.0 is provided

    @pytest.mark.skipif(
        not shutil.which("augustus"),
        reason="Augustus is not installed or not in PATH",
    )
    @pytest.mark.skipif(
        not shutil.which("miniprot"),
        reason="Miniprot is not installed or not in PATH",
    )
    def test_api_workflow(self, test_data_dir, temp_dir):
        """Test the complete buscolite API workflow."""
        # Set up paths
        genome_path = os.path.join(test_data_dir, "test_genome.fasta")
        proteins_path = os.path.join(test_data_dir, "test_proteins.fasta")
        lineage_path = os.path.join(test_data_dir, "mock_lineage")

        # Test genome mode
        genome_results, genome_missing, genome_stats, genome_config = runbusco(
            input=genome_path,
            lineage=lineage_path,
            mode="genome",
            species="anidulans",
            cpus=1,
            offset=2000,
            verbosity=3,
            check_augustus=False,  # Skip augustus check for testing
        )

        # Check genome mode results
        assert isinstance(genome_results, dict)
        assert isinstance(genome_missing, list)
        assert isinstance(genome_stats, dict)
        assert isinstance(genome_config, dict)

        # Test proteins mode
        proteins_results, proteins_missing, proteins_stats, proteins_config = runbusco(
            input=proteins_path,
            lineage=lineage_path,
            mode="proteins",
            species="anidulans",
            cpus=1,
            offset=2000,
            verbosity=3,
        )

        # Check proteins mode results
        assert isinstance(proteins_results, dict)
        assert isinstance(proteins_missing, list)
        assert isinstance(proteins_stats, dict)
        assert isinstance(proteins_config, dict)

        # Compare results
        assert genome_config == proteins_config
        assert genome_stats["total"] == proteins_stats["total"]
