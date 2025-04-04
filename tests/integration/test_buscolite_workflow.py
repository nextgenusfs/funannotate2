"""
Integration tests for the buscolite workflow.
"""

import os
import sys
import json
import tempfile
import shutil
import pytest
from pathlib import Path
from buscolite.busco import runbusco


class TestBUSCOliteWorkflow:
    """Test the buscolite workflow."""

    @pytest.fixture(scope="class")
    def test_data_dir(self):
        """Return the path to the test data directory."""
        return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

    @pytest.fixture(scope="function")
    def temp_dir(self):
        """Create a temporary directory for test outputs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield tmpdir

    @pytest.mark.skipif(
        not shutil.which("augustus"),
        reason="Augustus is not installed or not in PATH",
    )
    @pytest.mark.skipif(
        not shutil.which("miniprot"),
        reason="Miniprot is not installed or not in PATH",
    )
    def test_genome_mode(self, test_data_dir, temp_dir):
        """Test the buscolite workflow in genome mode."""
        # Set up paths
        genome_path = os.path.join(test_data_dir, "test_genome.fasta")
        lineage_path = os.path.join(test_data_dir, "mock_lineage")
        output_prefix = os.path.join(temp_dir, "test_output")

        # Run buscolite
        results, missing, stats, config = runbusco(
            input=genome_path,
            lineage=lineage_path,
            mode="genome",
            species="anidulans",
            cpus=1,
            offset=2000,
            verbosity=3,
            check_augustus=False,  # Skip augustus check for testing
        )

        # Check that the results are as expected
        assert isinstance(results, dict)
        assert isinstance(missing, list)
        assert isinstance(stats, dict)
        assert isinstance(config, dict)

        # Check that the stats are as expected
        assert "total" in stats
        assert "single-copy" in stats
        assert "fragmented" in stats
        assert "duplicated" in stats
        assert "missing" in stats

        # Check that the config is as expected
        assert config["name"] == "mock_lineage"
        assert config["species"] == "anidulans"
        assert config["domain"] == "eukaryota"

        # Write output files
        gff_path = f"{output_prefix}.buscolite.gff3"
        summary_path = f"{output_prefix}.buscolite.tsv"
        json_path = f"{output_prefix}.buscolite.json"

        # Write GFF file
        with open(gff_path, "w") as outfile:
            from buscolite.gff import gffwriter

            gffwriter(results, outfile)

        # Write summary file
        with open(summary_path, "w") as outfile:
            from buscolite.utilities import summary_writer

            summary_writer(results, missing, ["test"], config, outfile, mode="genome")

        # Write JSON file
        with open(json_path, "w") as outfile:
            outfile.write(json.dumps(results, indent=2))

        # Check that the output files exist
        assert os.path.exists(gff_path)
        assert os.path.exists(summary_path)
        assert os.path.exists(json_path)

        # Check that the GFF file is not empty
        assert os.path.getsize(gff_path) > 0

        # Check that the summary file is not empty
        assert os.path.getsize(summary_path) > 0

        # Check that the JSON file is not empty
        assert os.path.getsize(json_path) > 0

    @pytest.mark.skip(reason="buscolite is not installed")
    def test_proteins_mode(self, test_data_dir, temp_dir):
        """Test the buscolite workflow in proteins mode."""
        # Set up paths
        proteins_path = os.path.join(test_data_dir, "test_proteins.fasta")
        lineage_path = os.path.join(test_data_dir, "mock_lineage")
        output_prefix = os.path.join(temp_dir, "test_output")

        # Run buscolite
        results, missing, stats, config = runbusco(
            input=proteins_path,
            lineage=lineage_path,
            mode="proteins",
            species="anidulans",
            cpus=1,
            offset=2000,
            verbosity=3,
        )

        # Check that the results are as expected
        assert isinstance(results, dict)
        assert isinstance(missing, list)
        assert isinstance(stats, dict)
        assert isinstance(config, dict)

        # Check that the stats are as expected
        assert "total" in stats
        assert "single-copy" in stats
        assert "fragmented" in stats
        assert "duplicated" in stats
        assert "missing" in stats

        # Check that the config is as expected
        assert config["name"] == "mock_lineage"
        assert config["species"] == "anidulans"
        assert config["domain"] == "eukaryota"

        # Write output files
        summary_path = f"{output_prefix}.buscolite.tsv"
        json_path = f"{output_prefix}.buscolite.json"

        # Write summary file
        with open(summary_path, "w") as outfile:
            from buscolite.utilities import summary_writer

            summary_writer(results, missing, ["test"], config, outfile, mode="proteins")

        # Write JSON file
        with open(json_path, "w") as outfile:
            outfile.write(json.dumps(results, indent=2))

        # Check that the output files exist
        assert os.path.exists(summary_path)
        assert os.path.exists(json_path)

        # Check that the summary file is not empty
        assert os.path.getsize(summary_path) > 0

        # Check that the JSON file is not empty
        assert os.path.getsize(json_path) > 0
