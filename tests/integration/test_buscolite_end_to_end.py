"""
End-to-end integration tests for buscolite.
"""
import os
import sys
import json
import subprocess
import tempfile
import shutil
import pytest
from pathlib import Path


class TestBUSCOliteEndToEnd:
    """Test the buscolite end-to-end workflow."""

    @pytest.mark.skipif(
        not shutil.which("augustus"),
        reason="Augustus is not installed or not in PATH",
    )
    @pytest.mark.skipif(
        not shutil.which("miniprot"),
        reason="Miniprot is not installed or not in PATH",
    )
    def test_genome_mode_end_to_end(self, test_data_dir, temp_dir):
        """Test the buscolite end-to-end workflow in genome mode."""
        # Set up paths
        genome_path = os.path.join(test_data_dir, "test_genome.fasta")
        lineage_path = os.path.join(test_data_dir, "mock_lineage")
        output_prefix = os.path.join(temp_dir, "test_output")

        # Run buscolite
        result = subprocess.run(
            [
                "buscolite",
                "-i",
                genome_path,
                "-o",
                output_prefix,
                "-m",
                "genome",
                "-l",
                lineage_path,
                "-c",
                "1",
                "-s",
                "anidulans",
                "-f",
                "2000",
            ],
            capture_output=True,
            text=True,
        )

        # Check that the command ran successfully
        assert result.returncode == 0, f"Command failed with error: {result.stderr}"

        # Check that the output files exist
        gff_path = f"{output_prefix}.buscolite.gff3"
        summary_path = f"{output_prefix}.buscolite.tsv"
        json_path = f"{output_prefix}.buscolite.json"

        assert os.path.exists(gff_path), f"GFF file {gff_path} does not exist"
        assert os.path.exists(
            summary_path
        ), f"Summary file {summary_path} does not exist"
        assert os.path.exists(json_path), f"JSON file {json_path} does not exist"

        # Check that the output files are not empty
        assert os.path.getsize(gff_path) > 0, f"GFF file {gff_path} is empty"
        assert (
            os.path.getsize(summary_path) > 0
        ), f"Summary file {summary_path} is empty"
        assert os.path.getsize(json_path) > 0, f"JSON file {json_path} is empty"

        # Check the content of the JSON file
        with open(json_path, "r") as f:
            results = json.load(f)
            assert isinstance(results, dict)

        # Check the content of the summary file
        with open(summary_path, "r") as f:
            summary = f.read()
            assert "# BUSCO version:" in summary
            assert "# Lineage:" in summary
            assert "# Mode: genome" in summary

    def test_proteins_mode_end_to_end(self, test_data_dir, temp_dir):
        """Test the buscolite end-to-end workflow in proteins mode."""
        # Set up paths
        proteins_path = os.path.join(test_data_dir, "test_proteins.fasta")
        lineage_path = os.path.join(test_data_dir, "mock_lineage")
        output_prefix = os.path.join(temp_dir, "test_output")

        # Run buscolite
        result = subprocess.run(
            [
                "buscolite",
                "-i",
                proteins_path,
                "-o",
                output_prefix,
                "-m",
                "proteins",
                "-l",
                lineage_path,
                "-c",
                "1",
            ],
            capture_output=True,
            text=True,
        )

        # Check that the command ran successfully
        assert result.returncode == 0, f"Command failed with error: {result.stderr}"

        # Check that the output files exist
        summary_path = f"{output_prefix}.buscolite.tsv"
        json_path = f"{output_prefix}.buscolite.json"

        assert os.path.exists(
            summary_path
        ), f"Summary file {summary_path} does not exist"
        assert os.path.exists(json_path), f"JSON file {json_path} does not exist"

        # Check that the output files are not empty
        assert (
            os.path.getsize(summary_path) > 0
        ), f"Summary file {summary_path} is empty"
        assert os.path.getsize(json_path) > 0, f"JSON file {json_path} is empty"

        # Check the content of the JSON file
        with open(json_path, "r") as f:
            results = json.load(f)
            assert isinstance(results, dict)

        # Check the content of the summary file
        with open(summary_path, "r") as f:
            summary = f.read()
            assert "# BUSCO version:" in summary
            assert "# Lineage:" in summary
            assert "# Mode: proteins" in summary
